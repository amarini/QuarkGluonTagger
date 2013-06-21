#!/bin/bash

FILE=../data/SystDatabase.txt
mkdir -p ~/work/PlotQG/

do_dijet()
{
linenumber=0
cat ${FILE} | grep -v '^#' | while read varname ptmin ptmax rhomin rhomax etamin etamax alpha beta lmin lmax ; do
((linenumber++))
[ $etamax -eq 5 ] && etamax=4.7;
echo "Executing line DiJet $linenumber"
root  -l -b <<EOF
.L ComputeDoubleMinDiJet.C+
TCanvas *c;
Check(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${alpha},${beta},"${varname}",${lmin},${lmax},&c);
c->SaveAs("~/work/PlotQG/DiJetSystValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.pdf")
.q
EOF
done
}

do_zjet()
{
linenumber=0
cat ${FILE} | grep -v '^#' | while read varname ptmin ptmax rhomin rhomax etamin etamax alpha beta lmin lmax ; do
((linenumber++))
[ $etamax -eq 5 ] && etamax=4.7;
echo "Executing line ZJet $linenumber: varname=${varname} pt=[${ptmin},${ptmax}] rho=[${rhomin},${rhomax}] eta=[$etamin,$etamax] par=[$alpha,$beta,$lmin,$lmax]"
root  -l -b <<EOF
.L ComputeDoubleMin.C+
TCanvas *c;
Check(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${alpha},${beta},"${varname}",${lmin},${lmax},&c);
c->SaveAs("~/work/PlotQG/ZJetSystValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.pdf")
.q
EOF
done
}

do_zjet2()
{
linenumber=0
cat ${FILE} | grep -v '^#' | while read varname ptmin ptmax rhomin rhomax etamin etamax alpha beta lmin lmax ; do
((linenumber++))
[ $etamax -eq 5 ] && etamax=4.7;
echo "Executing line ZJet2 $linenumber"
root  -l -b <<EOF
.L ComputeDoubleMinZJet2.C+
TCanvas *c;
Check(${ptmin},${ptmax},${rhomin},${rhomax},${etamin},${etamax},${alpha},${beta},"${varname}",${lmin},${lmax},&c);
c->SaveAs("~/work/PlotQG/ZJetHSSystValidation_${varname}_pt${ptmin}_${ptmax}_eta${etamin}_${etamax}_rho${rhomin}_${rhomax}.pdf")
.q
EOF
done
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
