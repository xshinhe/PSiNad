./psinad -w -p param_FMO_NAFTW1_init5.json
./psinad -w -p param_FMO_NAFTW1_init5.json -load=default/record-dump0-100.ds:continue -d continue
./psinad -w -p param_FMO_NAFTW1_init5.json -load=default/record-dump0-100.ds:restart  -d restart
./psinad -w -p param_FMO_NAFTW1_init5_rst.json -load=./default/record-dump0-100.ds:restart -d restart2

