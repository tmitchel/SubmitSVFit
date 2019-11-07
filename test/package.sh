#!/bin/bash
pushd $CMSSW_BASE/src
tar  czf svfitter_packge.tar.gz HTT-utilities SubmitSVFit TauAnalysis 
xrdcp -f svfitter_packge.tar.gz root://cmseos.fnal.gov//store/user/$1/svfitter_package/code.tar.gz
rm svfitter_packge.tar.gz
popd


