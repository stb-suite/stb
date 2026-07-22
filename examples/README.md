# Examples

Hands-on, readable tutorials aimed at people **learning to use** stb_suite —
as opposed to `test/`, which is the developers' smoke-test suite (edge
cases, expected errors, PASS/FAIL counts). Each folder here is a small
lesson: what the tool does, why you'd reach for it, and a commented script
you run and read at the same time.

## How it's organized

Each example lives in a folder named after the same dotted code the
interactive `stb-suite` menu uses to jump straight to a tool (its main
prompt accepts codes like `1.1`, `2.3`, `4.1.2`). So `1.1-stb-inputfile/`
corresponds to typing `1.1` in the `stb-suite` menu, or running
`stb-inputfile` directly from the command line — both are shown in every
example. See `CLAUDE.md` (the `INPUT_TOOLS`, `STRUCTURE_TOOLS`,
`ANALYSIS_TOOLS`, `WORKFLOW_TOOLS`, `MLSIM_TOOLS`, `UTILITY_TOOLS`
dictionaries) for the full list of codes.

This is an ongoing effort — new examples are added a few at a time, one per
tool, not all at once.

## Two ways to run any tool

Every example shows both:

1. **Direct CLI** — call the `stb-*` command yourself with flags. Faster,
   scriptable, what you'd use once you know the tool.
2. **Interactive `stb-suite` menu** — run `stb-suite` and either navigate the
   category menus or type the tool's dotted code (e.g. `1.1`) straight from
   the main prompt. It walks you through the same choices as guided
   questions instead of flags — a good way to discover what a tool can do
   before memorizing its options.

Both paths call the exact same underlying tool and produce the exact same
output.

## Prerequisites

Install the package in editable mode before running any example:

```bash
cd stb-suite
pip install -e .
```

This puts every `stb-*` command (and `stb-suite` itself) on your PATH.

## Index

| Code | Tool            | Folder                                    |
|------|-----------------|--------------------------------------------|
| 1.1  | `stb-inputfile` | [`1.1-stb-inputfile/`](1.1-stb-inputfile/)  |
| 1.2  | `stb-kgrid`     | [`1.2-stb-kgrid/`](1.2-stb-kgrid/)          |
| 1.3  | `stb-kpath`     | [`1.3-stb-kpath/`](1.3-stb-kpath/)          |
| 1.4  | `stb-dftu`      | [`1.4-stb-dftu/`](1.4-stb-dftu/)            |
| 1.5  | `stb-fetch`     | [`1.5-stb-fetch/`](1.5-stb-fetch/)          |
| 1.6  | `stb-mlrelax`   | [`1.6-stb-mlrelax/`](1.6-stb-mlrelax/)      |
