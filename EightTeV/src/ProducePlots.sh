#!/bin/bash

FILE=../data/SystDatabase.txt
FILE=./SystZJet.txt
FILE=./SystDoubleMinZJetFineBins.txt
#FILE=./ProducePlot.txt
mkdir -p ~/work/PlotQG/


template() # parameters: $1=TYPE=DiJet,ZJet,ZJet2 
{
linenumber=0
type1=$1;
type2=$1;
[ $type2 == "ZJet" ] && type2="";
cat ${FILE} | grep -v '^#' | grep -v '^$' | while read varname ptmin ptmax rhomin rhomax etamin etamax alpha beta lmin lmax ; do
((linenumber++))
[ "$etamax" == "5" ] && etamax="4.7";
[ "$etamax" == "5.0" ] && etamax="4.7";
echo "Executing line $type1 $linenumber - ComputeDoubleMin$type2"
root  -l -b <<EOF
.L ComputeDoubleMin${type2}.C+
TCanvas *c;
Check(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${alpha},${beta},"${varname}",${lmin},${lmax},&c);
c->SaveAs("~/work/PlotQG/${type1}_SystValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.pdf")
c->SaveAs("~/work/PlotQG/${type1}_SystValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.root")
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

check_parsing()
{
cat ${FILE} | grep -v '^#' | while read varname ptmin ptmax rhomin rhomax etamin etamax alpha beta lmin lmax ; do
((linenumber++))
[ $etamax -eq 5 ] && etamax=4.7;
echo "Executing line Check $linenumber: varname=${varname} pt=[${ptmin},${ptmax}] rho=[${rhomin},${rhomax}] eta=[$etamin,$etamax] par=[$alpha,$beta,$lmin,$lmax]"
echo "Check(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${alpha},${beta},\"${varname}\",${lmin},${lmax},&c);"
done

}

check_parsing;
do_zjet;
#do_zjet2;
#do_dijet;
