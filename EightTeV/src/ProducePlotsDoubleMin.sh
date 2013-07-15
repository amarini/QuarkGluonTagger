#!/bin/bash

#FILE=./SystDoubleMinZJet.txt
#FILE=./SystDoubleMinZJetFineBins.txt
FILE=./SystDoubleZJetHbb.txt
#FILE=./ProducePlot.txt

DIR=~/work/PlotQG/Double3/
mkdir -p ${DIR}

template() # parameters: $1=TYPE=DiJet,ZJet,ZJet2,ZJetHbb
{
linenumber=0
type1=$1;
type2=$1;
[ $type2 == "ZJet" ] && type2="";
cat ${FILE} | grep -v '^#' | grep -v '^$' | while read varname ptmin ptmax rhomin rhomax etamin etamax a_q b_q a_g b_g lmin lmax ; do
((linenumber++))
[ "$etamax" == "5" ] && etamax="4.7";
[ "$etamax" == "5.0" ] && etamax="4.7";
echo "Executing line $type1 $linenumber - ComputeDoubleMin$type2 - [$ptmin,$ptmax][$rhomin,$rhomax][$etamin,$etamax][$a_q,$b_q,$a_g,$b_g,$lmin,$lmax]"
	[ "$type2" == "ZJet2" ] && [ $ptmax -lt 50 ] && continue;
root  -l -b <<EOF
.L ComputeDoubleMin${type2}.C+
TCanvas *c;
CheckDouble(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${a_q},${b_q},${a_g},${b_g},"${varname}",${lmin},${lmax},&c);
c->SaveAs("${DIR}/${type1}_SystDoubleValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.pdf")
c->SaveAs("${DIR}/${type1}_SystDoubleValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.root")
.q
EOF
done

}

do_dijet()
{
template DiJet ;
}

do_zjet()
{
template ZJet ;
}

do_zjet2()
{
template ZJet2 ;
}

do_zjetHbb()
{
template ZJetHbb ;
}

check_parsing()
{
cat ${FILE} | grep -v '^#' | while read varname ptmin ptmax rhomin rhomax etamin etamax a_q b_q a_g b_g lmin lmax ; do
((linenumber++))
[ "$etamax" == "5" ] && etamax=4.7;
[ "$etamax" == "5.0" ] && etamax=4.7;
echo "Executing line Check $linenumber: varname=${varname} pt=[${ptmin},${ptmax}] rho=[${rhomin},${rhomax}] eta=[$etamin,$etamax] par=[$a_q,$b_q,$a_g,$b_g,$lmin,$lmax]"
echo "CheckDouble(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${a_q},${b_q},${a_g},${b_g},\"${varname}\",${lmin},${lmax},&c);"
done

}

check_parsing;
do_zjetHbb;
do_zjet;
do_zjet2;
do_dijet;
