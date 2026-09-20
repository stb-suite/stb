# Bug (metodológico, não numérico): `stb-hirshfeldIons`/`stb-hirshfeldAnalysis` tentam construir uma referência de cátion para hidrogênio que é fisicamente impossível de calcular via SIESTA (H⁺ = próton nu, 0 elétrons)

**Encontrado em:** 19-20/09/2026, ao rodar o workflow Hirshfeld-I (item 20 do
menu, `stb-hirshfeldPrep`→`stb-hirshfeldIons`→`stb-hirshfeldAnalysis`) sobre as
16 combinações dopante×gás do projeto `calcogenides-GY-adsorption` (grafino
tipo-γ dopado com calcogênios + gás adsorvido). Os 4 sistemas com NH3 como
adsorvato (`O_NH3`, `Se_NH3`, `S_NH3`, `Te_NH3` — os únicos 4 com hidrogênio na
composição) falharam no Estágio 3 com:

```
[ERROR] Could not build a cation profile for species 'H_ads' from
'hirshfeld_charges/O_NH3/hirshfeld_study/ions/H_ads/cation': no '*.RHO' file
found in '.../ions/H_ads/cation' -- has SIESTA been run there yet (with
SaveRho true)? -- has SIESTA been run there yet?
```

Investigando a pasta `ions/H_ads/cation/`, o SIESTA nunca chegou a terminar:

```
new_DM -- step:     1
Initializing Density Matrix...
Note: For starting DM, Qtot, Tr[D*S] =          0.00000000          1.00000000
...
stepf: Fermi-Dirac step function
 Fermid: Iteration has not converged.
 Fermid: qtot,sumq=   0.0000000000000000        1.0372736004734628E-002
Fermid: Iteration has not converged.
Stopping Program from Node:    0
```

## Causa raiz

`stb-hirshfeldIons` (Stage 2) escreve incondicionalmente uma pasta
`ions/<espécie>/cation/` com `NetCharge +1.0` para **toda** espécie do sistema
combinado, sem checar `Z_val` (a carga de valência da espécie, já lida do
`.out` combinado ou da tabela de fallback — a mesma resolução usada pelo
`stb-bader`). Para hidrogênio, `Z_val = 1`: uma carga líquida de `+1.0`
significa **0 elétrons** — um próton nu em uma caixa de vácuo, não um átomo
ionizado no sentido usual.

Isso não é um corner-case numérico raro: é um problema **fisicamente mal
posto**. O solver de ocupação de Fermi-Dirac do SIESTA distribui elétrons
entre orbitais Kohn-Sham; com 0 elétrons não há nada para distribuir, e a
iteração de busca do nível de Fermi nunca converge (`qtot=0` no log acima), o
SIESTA aborta o programa. Nenhum ajuste de parâmetros de convergência (mistura,
`ElectronicTemperature`, etc.) resolve isso — o problema é estrutural, não de
tolerância numérica.

## Por que só hidrogênio é afetado

Hidrogênio é o único elemento comum com `Z_val = 1` nos pseudopotenciais
usuais (todos os outros elementos deste projeto — C, O, S, Se, Te, N — têm
`Z_val ≥ 4`, e mesmo elementos alcalinos raramente aparecem como espécie
adsorvida isolada num workflow de adsorção). Qualquer outra espécie com
`Z_val ≤ 1` sofreria do mesmo problema, mas isso é extremamente raro em
sistemas reais de SIESTA/DFT.

## Por que a referência "correta" não precisa de nenhum cálculo SIESTA

A densidade eletrônica de um próton nu (H⁺) é **exatamente zero em todo o
espaço** — isso não é uma aproximação nem um limite assintótico difícil de
alcançar numericamente, é um fato exato: sem elétrons, não há densidade
eletrônica para calcular. `stb-hirshfeldAnalysis` usa o perfil do cátion como
um dos dois extremos (`interpolate_reference_density`, em
`core/hirshfeld.py`) para interpolar a referência de cada átomo conforme sua
carga parcial:

```python
frac = min(abs(charge), cap)
neutral_i = np.interp(target_r, r_neutral, rho_neutral, left=0.0, right=0.0)
ion_i = np.interp(target_r, r_ion, rho_ion, left=0.0, right=0.0)
return (1.0 - frac) * neutral_i + frac * ion_i
```

Bastando fornecer `rho_ion = 0` (no mesmo grid radial da referência neutra) em
vez de tentar ler um `.RHO` de uma pasta SIESTA que nunca vai existir, essa
interpolação já produz o comportamento fisicamente correto: conforme a carga
do H se aproxima de +1 (cátion puro), a referência tende suavemente a zero.

## Correção proposta (patch em anexo, `hirshfeld_bare_proton_cation_fix.patch`)

Dois arquivos, aplicável via `git apply` a partir da raiz do repositório
(`stb-suite/src/stb/`):

1. **`hirshfeld_ions.py`** (Stage 2): antes de escrever a pasta do cátion,
   checa `valence_source[sym] <= 1.0`. Se verdadeiro, **pula** a geração/
   execução SIESTA daquela pasta (`cation_folder = None`), grava
   `"cation_zero_density": true` no manifesto (`hirshfeld_ions_manifest.json`)
   e imprime uma mensagem `[OK]` explicando o motivo (não é um erro
   silencioso). A pasta `ions/<espécie>/anion/` continua sendo gerada
   normalmente — H⁻ (2 elétrons) é um sistema perfeitamente calculável.

2. **`hirshfeld_analysis.py`** (Stage 3): ao carregar `cation_profiles[sym]`,
   checa a flag `cation_zero_density` do manifesto antes de tentar ler
   `*.RHO`. Se marcada, reaproveita o grid radial (`r`) já lido do perfil
   neutro da mesma espécie e usa `rho = np.zeros_like(r)` — nenhuma leitura de
   arquivo, nenhuma dependência de um cálculo SIESTA que nunca existiu.

Ambos os arquivos sobem a versão de `2.0.0` para `2.1.0` (convenção já usada
neste módulo para mudanças de comportamento).

### Validação (projeto `calcogenides-GY-adsorption`, 4 sistemas com NH3)

Antes do patch: Stage 3 falhava com o erro acima nos 4 sistemas
(`O_NH3`, `Se_NH3`, `S_NH3`, `Te_NH3`) — 0/4 relatórios gerados.

Depois do patch (mesmos dados de entrada, nenhum cálculo SIESTA extra além do
ânion de H, que já era necessário): os 4 convergiram normalmente em 14-15
iterações (tolerância padrão 0.005 e⁻), carga total do sistema combinado
≈0 em todos (checagem de sanidade da própria ferramenta):

| Sistema | Iterações | Convergiu | Carga total (e⁻) | H_ads médio (e⁻) | N_ads (e⁻) |
|---|---:|---|---:|---:|---:|
| O_NH3  | 14 | sim | −0.0000 | +0.2745 | −0.8823 |
| Se_NH3 | 14 | sim | −0.0001 | +0.3069 | −0.9374 |
| S_NH3  | 14 | sim | −0.0001 | +0.2889 | −0.9139 |
| Te_NH3 | 15 | sim | −0.0001 | +0.3304 | −0.9371 |

Os valores de H (levemente positivo, ~+0.27 a +0.33 e⁻, "cation-like" nos 4
sistemas) e N (fortemente negativo, ~−0.88 a −0.94 e⁻) são quimicamente
sensatos para N-H em uma amina — o N mais eletronegativo puxa densidade dos 3
H — e a carga líquida da molécula (N + 3H) fica próxima de zero/levemente
negativa em todos os 4 casos, consistente com o papel de doador de elétrons já
visto para NH3 via Δρ (seção 10 do relatório do projeto) nos outros 3
adsorbatos.

## Por que a correção é segura (não é um ajuste arbitrário)

- Não introduz nenhum parâmetro livre nem aproximação numérica — `ρ(r)=0`
  para H⁺ é um resultado exato, não uma escolha de conveniência.
- Não muda o comportamento de nenhuma outra espécie (`valence_source[sym] > 1`
  em todos os outros casos comuns) nem do próprio H⁻ (ânion, que continua
  sendo um cálculo SIESTA real de 2 elétrons).
- A extrapolação `left=0.0, right=0.0` que `interpolate_reference_density` já
  usa para qualquer perfil fora do seu próprio `r_max` significa que um
  `rho_ion` totalmente zero se comporta exatamente como o "vazio" que a
  função já trata como caso de borda — não há caminho de código novo sendo
  exercitado, só um valor de entrada diferente.
- É estritamente mais informativo que o comportamento anterior: antes, o
  workflow simplesmente abortava sem terminar nenhuma espécie com H;
  agora processa normalmente e ainda avisa explicitamente (`[OK] ... cation
  SKIPPED ...`) que aquela espécie usou a referência analítica, não uma
  omissão silenciosa.

## Arquivos deste diretório

- `hirshfeld_bare_proton_cation_fix.patch` — patch unificado (`git diff`
  contra o `main` do repositório `stb-suite/stb`, commit `8e65271`),
  aplicável com `git apply hirshfeld_bare_proton_cation_fix.patch` a partir da
  raiz do repositório.
