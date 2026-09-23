# Pathway and enrichment analysis

Transparent scaffold for pathway scoring and enrichment from gene/protein statistics. The portable runner validates identifiers, declares a tested universe, ranks pathways, and records audit metadata. The original GO/KEGG/Hallmark/Reactome R workflow remains under `src/legacy/`.

## Quick start

```bash
python workflow.py --input examples/input.csv --config config/workflow.json --output output
python -m unittest discover -s tests -v
```

For production ORA/GSEA, pin gene-set and annotation releases, document identifier conversion and universe selection, correct for multiple testing, and inspect coverage/leading features instead of interpreting pathway names alone.
