#!/bin/bash

# --- Configuração ---
# Nome do script Python que queremos testar
SCRIPT_NAME="stb-translate"
# Comando para executar o Python (pode ser "python" ou "python3")
PYTHON_CMD=""
# Diretório para guardar os ficheiros de teste
TEST_DIR="test_files"

# Cores para o output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # Sem Cor

# --- Função de Verificação ---
# Esta função verifica se o comando anterior (a conversão) foi bem-sucedido
# e se o ficheiro de saída esperado foi criado e não está vazio.
check_success() {
    # $? é o código de saída do último comando. 0 = sucesso.
    # -s "$1" verifica se o ficheiro $1 existe E não está vazio.
    if [ $? -eq 0 ] && [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (Ficheiro '$1' criado)"
    else
        echo -e " ... ${RED}FALHA${NC}"
    fi
}

# --- 1. Preparação ---
echo "--- Iniciando Testador para STB-Translate ---"

# Cria o diretório de teste (e não dá erro se já existir)
mkdir -p $TEST_DIR
echo "Diretório de teste '$TEST_DIR' preparado."

# --- 2. Criação de Ficheiros de Exemplo ---
# Vamos criar uma estrutura simples de 2 átomos de Silício (Si)
# em vários formatos de ficheiro.

echo "Gerando ficheiros de entrada de exemplo..."

# --- Exemplo: si.cif ---
cat > $TEST_DIR/si.cif << 'EOF'
data_silicon
_cell_length_a     5.43
_cell_length_b     5.43
_cell_length_c     5.43
_cell_angle_alpha  90.0
_cell_angle_beta   90.0
_cell_angle_gamma  90.0
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.00 0.00 0.00
Si2 Si 0.25 0.25 0.25
EOF

# --- Exemplo: si.poscar ---
cat > $TEST_DIR/si.poscar << 'EOF'
Silicon Unit Cell
1.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
EOF

# --- Exemplo: si.fdf (O NOVO FORMATO DE ENTRADA) ---
cat > $TEST_DIR/si.fdf << 'EOF'
SystemName          Silicon
LatticeConstant     5.43 Ang

%block LatticeVectors
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
%endblock LatticeVectors

AtomicCoordinatesFormat  Fractional

%block ChemicalSpeciesLabel
1  14  Si
%endblock ChemicalSpeciesLabel

%block AtomicCoordinatesAndAtomicSpecies
0.00 0.00 0.00 1
0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

# --- Exemplo: si.siesta (Formato antigo) ---
cat > $TEST_DIR/si.siesta << 'EOF'
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
4:
1  14  0.00 0.00 0.00
1  14  0.25 0.25 0.25
EOF

# --- Exemplo: si.xyz ---
cat > $TEST_DIR/si.xyz << 'EOF'
2
Silicon Cartesian Coords
Si 0.0000 0.0000 0.0000
Si 1.3575 1.3575 1.3575
EOF

# --- Exemplo: si.dftb ---
cat > $TEST_DIR/si.dftb << 'EOF'
2 F
Si
1 1 0.00 0.00 0.00
2 1 0.25 0.25 0.25
0.0 0.0 0.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
EOF

# --- Exemplo: si.xsf ---
cat > $TEST_DIR/si.xsf << 'EOF'
CRYSTAL
PRIMVEC
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
PRIMCOORD
2 1
14 0.0000 0.0000 0.0000
14 1.3575 1.3575 1.3575
EOF

# --- Exemplo: si.fhi ---
cat > $TEST_DIR/si.fhi << 'EOF'
lattice_vector 5.430 0.000 0.000
lattice_vector 0.000 5.430 0.000
lattice_vector 0.000 0.000 5.430
atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
EOF

echo -e "${GREEN}Ficheiros de exemplo gerados em '$TEST_DIR/'${NC}"


# --- 3. Execução dos Testes ---

# Usamos --no-intro para um output mais limpo nos testes
BASE_CMD="$PYTHON_CMD $SCRIPT_NAME --no-intro"

# --- Teste de LEITORES (Todos os formatos -> CIF) ---
echo -e "\n--- Testando Leitores (Todos os formatos -> .cif) ---"

echo -n "Testando: POSCAR -> CIF"
$BASE_CMD -if poscar -i $TEST_DIR/si.poscar -of cif -o $TEST_DIR/out_from_poscar.cif
check_success "$TEST_DIR/out_from_poscar.cif"

echo -n "Testando: FDF -> CIF ${YELLOW}(NOVO!)${NC}"
$BASE_CMD -if fdf -i $TEST_DIR/si.fdf -of cif -o $TEST_DIR/out_from_fdf.cif
check_success "$TEST_DIR/out_from_fdf.cif"

echo -n "Testando: SIESTA -> CIF"
$BASE_CMD -if siesta -i $TEST_DIR/si.siesta -of cif -o $TEST_DIR/out_from_siesta.cif
check_success "$TEST_DIR/out_from_siesta.cif"

echo -n "Testando: XYZ (+lattice-vectors) -> CIF"
$BASE_CMD -if xyz -i $TEST_DIR/si.xyz --lattice-vectors 5.43 0 0 0 5.43 0 0 0 5.43 -of cif -o $TEST_DIR/out_from_xyz.cif
check_success "$TEST_DIR/out_from_xyz.cif"

echo -n "Testando: DFTB -> CIF"
$BASE_CMD -if dftb -i $TEST_DIR/si.dftb -of cif -o $TEST_DIR/out_from_dftb.cif
check_success "$TEST_DIR/out_from_dftb.cif"

echo -n "Testando: XSF -> CIF"
$BASE_CMD -if xsf -i $TEST_DIR/si.xsf -of cif -o $TEST_DIR/out_from_xsf.cif
check_success "$TEST_DIR/out_from_xsf.cif"

echo -n "Testando: FHI -> CIF"
$BASE_CMD -if fhi -i $TEST_DIR/si.fhi -of cif -o $TEST_DIR/out_from_fhi.cif
check_success "$TEST_DIR/out_from_fhi.cif"


# --- Teste de ESCRITORES (CIF -> Todos os formatos) ---
echo -e "\n--- Testando Escritores (.cif -> Todos os formatos) ---"

echo -n "Testando: CIF -> POSCAR"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of poscar -o $TEST_DIR/out_from_cif.poscar
check_success "$TEST_DIR/out_from_cif.poscar"

echo -n "Testando: CIF -> FDF"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of fdf -o $TEST_DIR/out_from_cif.fdf
check_success "$TEST_DIR/out_from_cif.fdf"

echo -n "Testando: CIF -> XYZ"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of xyz -o $TEST_DIR/out_from_cif.xyz
check_success "$TEST_DIR/out_from_cif.xyz"

echo -n "Testando: CIF -> DFTB"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of dftb -o $TEST_DIR/out_from_cif.dftb
check_success "$TEST_DIR/out_from_cif.dftb"

echo -n "Testando: CIF -> XSF"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of xsf -o $TEST_DIR/out_from_cif.xsf
check_success "$TEST_DIR/out_from_cif.xsf"

echo -n "Testando: CIF -> FHI"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of fhi -o $TEST_DIR/out_from_cif.fhi
check_success "$TEST_DIR/out_from_cif.fhi"

echo -n "Testando: CIF -> LAMMPS ${YELLOW}(NOVO!)${NC}"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of lammps -o $TEST_DIR/out_from_cif.data
check_success "$TEST_DIR/out_from_cif.data"

# Confirma que o formato lammps tem uma seção "Masses" (marca de um ficheiro
# lammps-data válido) e que a checagem de integridade (ida e volta) não falha
# -- essa reconferência usa ase.io.read(format='lammps-data') diretamente
# (lammps é só formato de SAÍDA, sem getatomsandvectors_lammps de entrada).
grep -q "Masses" $TEST_DIR/out_from_cif.data
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: Seção 'Masses' encontrada.${NC}"
else
    echo -e "   -> ${RED}Falha: Seção 'Masses' NÃO encontrada.${NC}"
fi


# --- Teste da Flag -cf (Conversão de Coordenadas) ---
echo -e "\n--- Testando Conversão de Coordenadas (-cf) ${YELLOW}(NOVO!)${NC} ---"

# O nosso si.poscar tem coordenadas 'Direct' (Fracionárias)
echo -n "Testando: POSCAR (Direct) -> FDF (Cartesian)"
$BASE_CMD -if poscar -i $TEST_DIR/si.poscar -of fdf -o $TEST_DIR/out_cart.fdf -cf cartesian
check_success "$TEST_DIR/out_cart.fdf"

# Vamos verificar se o ficheiro de saída contém 'Ang' (indicando Cartesiano)
grep -q "Ang" $TEST_DIR/out_cart.fdf
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: Formato 'Ang' (Cartesiano) encontrado.${NC}"
else
    echo -e "   -> ${RED}Falha: Formato 'Ang' (Cartesiano) não encontrado.${NC}"
fi

# O nosso si.cif é lido como Cartesiano
echo -n "Testando: CIF (Cartesian) -> POSCAR (Direct)"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of poscar -o $TEST_DIR/out_direct.poscar -cf direct
check_success "$TEST_DIR/out_direct.poscar"

# Vamos verificar se o ficheiro de saída contém 'Direct'
grep -q "Direct" $TEST_DIR/out_direct.poscar
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: Formato 'Direct' encontrado.${NC}"
else
    echo -e "   -> ${RED}Falha: Formato 'Direct' não encontrado.${NC}"
fi

# POSCAR com a linha de coordenadas em minúsculas ('direct'): VASP só olha a primeira letra,
# e o resultado sem -cf deve ser 'Fractional', não 'Ang' com valores fracionários.
echo -n "Testando: POSCAR ('direct' minúsculo) -> FDF sem -cf"
sed 's/^Direct$/direct/' $TEST_DIR/si.poscar > $TEST_DIR/si_lower.poscar
$BASE_CMD -if poscar -i $TEST_DIR/si_lower.poscar -of fdf -o $TEST_DIR/out_lower.fdf
check_success "$TEST_DIR/out_lower.fdf"
grep -q "AtomicCoordinatesFormat  Fractional" $TEST_DIR/out_lower.fdf
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: 'Fractional' preservado.${NC}"
else
    echo -e "   -> ${RED}Falha: coordenadas fracionárias não rotuladas como 'Fractional'.${NC}"
fi


# --- 3b. Formatos automáticos, relatório, verificação e avisos ---
echo -e "\n--- Testando formato automático (entrada e saída), relatório e verificações ---"
NPASS=0; NFAIL=0
ok()  { echo -e "   -> ${GREEN}Verificado:${NC} $1"; NPASS=$((NPASS+1)); }
bad() { echo -e "   -> ${RED}Falha:${NC} $1"; NFAIL=$((NFAIL+1)); }
has() { if grep -q -- "$1" "$2" 2>/dev/null; then ok "'$1' em $2"; else bad "'$1' NÃO está em $2"; fi; }

echo "Testando: sem -if/-of, os dois formatos saem dos nomes"
$BASE_CMD -i $TEST_DIR/si.poscar -o $TEST_DIR/auto_out.fdf --save-report > $TEST_DIR/log_auto.txt 2>&1
[ -s $TEST_DIR/auto_out.fdf ] && ok "auto_out.fdf escrito" || bad "auto_out.fdf não escrito"
has "poscar (auto, from the name 'si.poscar')" $TEST_DIR/log_auto.txt
has "fdf (auto, from the name 'auto_out.fdf')" $TEST_DIR/log_auto.txt
for sct in "RUN METADATA" "INPUT" "STRUCTURE" "SYMMETRY" "OUTPUT" "SUMMARY & FILES" "WARNINGS"; do
    has "$sct" $TEST_DIR/log_auto.txt
done
has "Verification (the written file read back)" $TEST_DIR/stb_translate_report.txt
for ext in xyz cif vasp gen in xsf data; do
    $BASE_CMD -i $TEST_DIR/si.poscar -o $TEST_DIR/auto_ext.$ext > $TEST_DIR/log_ext_$ext.txt 2>&1
    [ -s $TEST_DIR/auto_ext.$ext ] && ok "saída .$ext detectada e escrita" || bad "saída .$ext não escrita"
done
echo "Testando: extensão desconhecida pede -of"
$BASE_CMD -i $TEST_DIR/si.poscar -o $TEST_DIR/auto_ext.unknown > $TEST_DIR/log_unknown.txt 2>&1
[ $? -eq 2 ] && ok "erro limpo (código 2)" || bad "extensão desconhecida não deu erro"
has "Could not infer the output format" $TEST_DIR/log_unknown.txt

echo "Testando: ida e volta .gen numa célula não cúbica (a e c eram trocados na leitura)"
python3 - <<PY
from ase.build import hcp0001
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
s = hcp0001("Mg", size=(2, 2, 3), vacuum=6.0)
structure_io.write_fdf(structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(s), coord_format="fractional"),
                       "$TEST_DIR/mg_slab.fdf")
PY
$BASE_CMD -i $TEST_DIR/mg_slab.fdf -o $TEST_DIR/mg_slab.gen > $TEST_DIR/log_gen1.txt 2>&1
$BASE_CMD -i $TEST_DIR/mg_slab.gen -o $TEST_DIR/mg_back.fdf > $TEST_DIR/log_gen2.txt 2>&1
has "largest atom displacement (Ang) |         | 0.00e+00 | ok" $TEST_DIR/log_gen1.txt
has "largest atom displacement (Ang) |         | 0.00e+00 | ok" $TEST_DIR/log_gen2.txt
python3 - <<PY
import sys, numpy as np
from stb import translate as t
a = t.build_ase_atoms(*t.read_structure("fdf", "$TEST_DIR/mg_slab.fdf"))
b = t.build_ase_atoms(*t.read_structure("fdf", "$TEST_DIR/mg_back.fdf"))
sys.exit(0 if np.allclose(a.cell, b.cell, atol=1e-6) and np.allclose(a.positions, b.positions, atol=1e-5) else 1)
PY
[ $? -eq 0 ] && ok "fdf -> gen -> fdf devolve a mesma célula e posições" || bad "a ida e volta .gen mudou a estrutura"

echo "Testando: coordenadas cartesianas rotuladas como Direct geram aviso"
python3 - <<PY
from ase.io import read
a = read("$TEST_DIR/si.cif")
lines = ["bad", "1.0"] + [" ".join(f"{x:.6f}" for x in v) for v in a.cell] + ["Si", str(len(a)), "Direct"] + \
        [" ".join(f"{x:.6f}" for x in p) for p in a.positions]
open("$TEST_DIR/bad_mixup.vasp", "w").write("\n".join(lines) + "\n")
PY
$BASE_CMD -i $TEST_DIR/bad_mixup.vasp --dry-run > $TEST_DIR/log_mixup.txt 2>&1
has "a Cartesian file labelled Direct/Fractional?" $TEST_DIR/log_mixup.txt

echo "Testando: CIF simetrizado (padrão) vs --cif-symprec 0 (P1, como lido)"
python3 - <<PY
import numpy as np
from ase.build import bulk
from pymatgen.io.ase import AseAtomsAdaptor
from stb.core import structure_io
s = bulk("Si", "diamond", a=5.43, cubic=True)
s.positions += np.random.default_rng(7).uniform(-0.002, 0.002, s.positions.shape)
structure_io.write_fdf(structure_io.from_pymatgen(AseAtomsAdaptor.get_structure(s), coord_format="fractional"),
                       "$TEST_DIR/si_noisy.fdf")
PY
$BASE_CMD -i $TEST_DIR/si_noisy.fdf -o $TEST_DIR/si_sym.cif > $TEST_DIR/log_cifsym.txt 2>&1
$BASE_CMD -i $TEST_DIR/si_noisy.fdf -o $TEST_DIR/si_p1.cif --cif-symprec 0 > $TEST_DIR/log_cifp1.txt 2>&1
has "| symmetrized" $TEST_DIR/log_cifsym.txt
has "_symmetry_space_group_name_H-M   'P 1'" $TEST_DIR/si_p1.cif
has "largest atom displacement (Ang) |         | 0.00e+00 | ok" $TEST_DIR/log_cifp1.txt

echo "Testando: POSCAR com fator de escala negativo (= volume da célula)"
python3 - <<PY
from ase.io import read
a = read("$TEST_DIR/si.cif")
v = a.cell.volume
lines = ["neg", f"{-v:.6f}"] + [" ".join(f"{x / v ** (1/3):.10f}" for x in c) for c in a.cell] + \
        ["Si", str(len(a)), "Direct"] + [" ".join(f"{x:.8f}" for x in p) for p in a.get_scaled_positions()]
open("$TEST_DIR/si_negscale.vasp", "w").write("\n".join(lines) + "\n")
PY
VOL=$(python3 -c "from ase.io import read; print(f'{read(\"$TEST_DIR/si.cif\").cell.volume:.4f}')")
$BASE_CMD -i $TEST_DIR/si_negscale.vasp --dry-run > $TEST_DIR/log_negscale.txt 2>&1
has "Volume            : $VOL Ang^3" $TEST_DIR/log_negscale.txt

echo "Testando: lote com formatos de entrada diferentes, um arquivo inválido"
mkdir -p $TEST_DIR/batch_in && cp $TEST_DIR/si.poscar $TEST_DIR/si.cif $TEST_DIR/si.fdf $TEST_DIR/batch_in/ \
    && printf 'nada\n' > $TEST_DIR/batch_in/broken.cif
$BASE_CMD -i "$TEST_DIR/batch_in/*" -of xyz --out-dir $TEST_DIR/batch_out > $TEST_DIR/log_batch.txt 2>&1
[ $? -eq 1 ] && ok "lote termina com código 1 (um arquivo falhou)" || bad "código de saída do lote"
has "1 of 4 FAILED" $TEST_DIR/log_batch.txt
[ -s $TEST_DIR/batch_out/si.xyz ] && ok "os arquivos bons foram escritos" || bad "lote não escreveu os bons"

echo "Testando: menu (stb-suite 3.1), formato de saída automático (opção 0)"
printf '3.1\n%s\n0\n%s\n0\n3\nn\nn\n\n0\n' "$TEST_DIR/si.poscar" "$TEST_DIR/menu_out.cif" \
    | stb-suite > $TEST_DIR/log_menu.txt 2>&1
[ -s $TEST_DIR/menu_out.cif ] && ok "menu escreveu menu_out.cif" || bad "menu não escreveu a saída"
has "Output format:" $TEST_DIR/log_menu.txt

echo -e "Novos testes: ${GREEN}$NPASS ok${NC}, ${RED}$NFAIL falhas${NC}"


# --- 4. Limpeza ---
echo -e "\n--- Testes Concluídos ---"
read -p "Deseja remover o diretório '$TEST_DIR' e todos os ficheiros de teste? (s/n) " -n 1 -r
echo # Mover para a próxima linha
if [[ $REPLY =~ ^[Ss]$ ]]; then
    echo "A limpar ficheiros de teste..."
    rm -rf $TEST_DIR
    echo -e "${GREEN}Limpeza concluída.${NC}"
else
    echo "Os ficheiros de teste foram mantidos em '$TEST_DIR/' para inspeção."
fi
