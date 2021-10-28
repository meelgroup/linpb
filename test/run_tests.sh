#!/bin/bash

time=$1
logfolder="/tmp/roundingsat/$2"
binary=$3
options=$4


SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

errors=0
tested=0
completed=0

echo "###########################"
echo "########## START ##########"
echo "###########################"
echo ""
echo "timeout: $time"
echo "data: $logfolder"
echo "binary: $binary"
echo "options: $options"
echo ""

declare -a arr_dec=(
"cnf/ec-rand4regsplit-v030-n1.cnf*UNSATISFIABLE"
"opb/dec/256ebits_any.opb*UNSATISFIABLE"
"opb/dec/256ebits_0.opb*UNSATISFIABLE"
"opb/dec/128ebits_any.opb*UNSATISFIABLE"
"opb/dec/128ebits_0.opb*UNSATISFIABLE"
"opb/dec/22array_alg_ineq7.opb*UNSATISFIABLE"
"opb/dec/21array_alg_ineq7.opb*UNSATISFIABLE"
"opb/dec/32array_alg_ineq5.opb*UNSATISFIABLE"
"opb/dec/air01.0.s.opb*SATISFIABLE"
"opb/dec/air01.0.u.opb*UNSATISFIABLE"
"opb/dec/air02.0.s.opb*SATISFIABLE"
"opb/dec/air02.0.u.opb*UNSATISFIABLE"
"opb/dec/air06.0.s.opb*SATISFIABLE"
"opb/dec/air06.0.u.opb*UNSATISFIABLE"
"opb/dec/bm23.0.s.opb*SATISFIABLE"
"opb/dec/bm23.0.u.opb*UNSATISFIABLE"
"opb/dec/cracpb1.0.s.opb*SATISFIABLE"
"opb/dec/cracpb1.0.u.opb*UNSATISFIABLE"
"opb/dec/diamond.0.d.opb*UNSATISFIABLE"
"opb/dec/lp4l.0.s.opb*SATISFIABLE"
"opb/dec/lp4l.0.u.opb*UNSATISFIABLE"
"opb/dec/p0040.0.s.opb*SATISFIABLE"
"opb/dec/p0040.0.u.opb*UNSATISFIABLE"
"opb/dec/p0291.0.s.opb*SATISFIABLE"
"opb/dec/p0291.0.u.opb*UNSATISFIABLE"
"opb/dec/pipex.0.s.opb*SATISFIABLE"
"opb/dec/pipex.0.u.opb*UNSATISFIABLE"
"opb/dec/sentoy.0.s.opb*SATISFIABLE"
"opb/dec/sentoy.0.u.opb*UNSATISFIABLE"
"opb/dec/stein9.0.s.opb*SATISFIABLE"
"opb/dec/stein9.0.u.opb*UNSATISFIABLE"
"opb/dec/stein15.0.s.opb*SATISFIABLE"
"opb/dec/stein15.0.u.opb*UNSATISFIABLE"
)

echo "########## simple #########"
echo ""

for j in "${arr_dec[@]}"; do
    formula="$(cut -d'*' -f1 <<<$j)"
    logfile="$logfolder/$formula"
    mkdir -p `dirname $logfile`
    echo -n "" > $logfile.proof
    echo -n "" > $logfile.formula
    formula="$SCRIPTPATH/instances/$formula"
    if [ ! -f "$formula" ]; then
        echo "$formula does not exist."
        exit 1
    fi
    obj="$(cut -d'*' -f2 <<<$j)"
    echo "running $binary $formula $options --proof-log=$logfile"
    output=`timeout $time $binary $formula $options --proof-log=$logfile 2>&1 | awk '/^o|SATISFIABLE|.*Assertion.*/ {print $2}'`
    if [ "$output" != "" ] && [ "$output" != "$obj" ]; then
        errors=`expr 1000 + $errors`
        echo "wrong output: $output vs $obj"
    fi
    echo "verifying $logfile"
    wc -l $logfile.proof
    veripb $logfile.formula $logfile.proof -d --arbitraryPrecision
    errors=`expr $? + $errors`
    echo $errors
    tested=`expr 1 + $tested`
    echo $tested
    echo ""
done

echo "########## arbitrary precision #########"
echo ""

for j in "${arr_dec[@]}"; do
    formula="$(cut -d'*' -f1 <<<$j)"
    logfile="$logfolder/$formula"
    mkdir -p `dirname $logfile`
    echo -n "" > $logfile.proof
    echo -n "" > $logfile.formula
    formula="$SCRIPTPATH/instances/$formula"
    if [ ! -f "$formula" ]; then
        echo "$formula does not exist."
        exit 1
    fi
    obj="$(cut -d'*' -f2 <<<$j)"
    echo "running $binary $formula $options --bits-learned=0 --bits-overflow=0 --bits-reduced=0 --proof-log=$logfile"
    output=`timeout $time $binary $formula $options --bits-learned=0 --bits-overflow=0 --bits-reduced=0 --proof-log=$logfile 2>&1 | awk '/^o|Error:|SATISFIABLE|.*Assertion.*/ {print $2}'`
    if [ "$output" != "" ] && [ "$output" != "$obj" ]; then
        errors=`expr 1000 + $errors`
        echo "wrong output: $output vs $obj"
    fi
    echo "verifying $logfile"
    wc -l $logfile.proof
    veripb $logfile.formula $logfile.proof -d --arbitraryPrecision
    errors=`expr $? + $errors`
    echo $errors
    tested=`expr 1 + $tested`
    echo $tested
    echo ""
done

declare -a arr_modes=(
linear
core-guided
core-guided
core-guided
hybrid
hybrid
hybrid
)

declare -a arr_lazy=(
simple
simple
lazy
reified
simple
lazy
reified
)

declare -a arr_opt=(
"wcnf/driverlog01bc.wcsp.dir.wcnf*2245"
"opb/opt/normalized-single-obj-f47-DC-Side1.seq-B-2-1-EDCBAir.opb*-1593213266"
"opb/opt/enigma.opb*0"
"opb/opt/stein9.opb*5"
"opb/opt/stein15.opb*9"
"opb/opt/stein27.opb*18"
"opb/opt/stein45.opb*30"
"opb/opt/p0033.opb*3089"
"opb/opt/p0040.opb*62027"
"opb/opt/p0201.opb*7615"
"opb/opt/p0282.opb*258411"
"opb/opt/p0291.opb*7609041"
"opb/opt/p0548.opb*8691"
"opb/opt/mod008.opb*307"
"opb/opt/mod010.opb*6548"
"opb/opt/air01.opb*6796"
"opb/opt/air02.opb*7810"
"opb/opt/air03.opb*340160"
"opb/opt/air06.opb*49649"
"opb/opt/pipex.opb*788263"
"opb/opt/sentoy.opb*-7772"
"opb/opt/bm23.opb*34"
"opb/opt/l152lav.opb*4722"
"opb/opt/lp4l.opb*2967"
"opb/opt/lseu.opb*1120"
"opb/opt/cracpb1.opb*22199"
)

for idx in "${!arr_modes[@]}"; do
    mode=${arr_modes[$idx]}
    lazy=${arr_lazy[$idx]}
    echo "########## $mode lazy=$lazy ##########"
    echo ""
    for j in "${arr_opt[@]}"; do
        formula="$(cut -d'*' -f1 <<<$j)"
        logfile="$logfolder/$formula"
        mkdir -p `dirname $logfile`
        echo -n "" > $logfile.proof
        echo -n "" > $logfile.formula
        formula="$SCRIPTPATH/instances/$formula"
        if [ ! -f "$formula" ]; then
            echo "$formula does not exist."
            exit 1
        fi
        obj="$(cut -d'*' -f2 <<<$j)"
        echo "running $binary $formula $options --opt-mode=$mode --cg-encoding=$lazy --proof-log=$logfile"
        output=`timeout $time $binary $formula $options --opt-mode=$mode --cg-encoding=$lazy --proof-log=$logfile 2>&1 | awk '/^o|Error:|UNSATISFIABLE|.*Assertion.*/ {print $2}'`
        if [ "$output" != "" ] && [ "$output" != "$obj" ]; then
            errors=`expr 1000 + $errors`
            echo "wrong output: $output vs $obj"
        fi
        echo "verifying $logfile"
        wc -l $logfile.proof
        veripb $logfile.formula $logfile.proof -d --arbitraryPrecision
        errors=`expr $? + $errors`
        echo $errors
        tested=`expr 1 + $tested`
        echo $tested
        echo ""
    done
done

echo "########## no proofs ##########"
echo ""
for j in "${arr_opt[@]}"; do
    formula="$(cut -d'*' -f1 <<<$j)"
    logfile="$logfolder/$formula"
    mkdir -p `dirname $logfile`
    echo -n "" > $logfile.proof
    echo -n "" > $logfile.formula
    formula="$SCRIPTPATH/instances/$formula"
    if [ ! -f "$formula" ]; then
        echo "$formula does not exist."
        exit 1
    fi
    obj="$(cut -d'*' -f2 <<<$j)"
    echo "running $binary $formula $options --opt-mode=hybrid"
    output=`timeout $time $binary $formula $options --opt-mode=hybrid 2>&1 | awk '/^o|Error:|UNSATISFIABLE|.*Assertion.*/ {print $2}'`
    if [ "$output" != "" ] && [ "$output" != "$obj" ]; then
        errors=`expr 1000 + $errors`
        echo "wrong output: $output vs $obj"
    fi
    echo $errors
    tested=`expr 1 + $tested`
    echo $tested
    echo ""
done

echo "tested: $tested"
echo "errors: $errors"

# command to remove soplex asserts:
# grep -rl "\sassert(.*" . | xargs sed -i 's/assert(/assert(true || /g'
