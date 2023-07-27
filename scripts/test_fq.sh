cairo-compile ./tests/test_fq.cairo --output test_fq_compiled.json

cairo-run --program=test_fq_compiled.json --print_output --layout=all