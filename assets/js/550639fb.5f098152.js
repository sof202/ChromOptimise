"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[6366],{8095:(e,i,t)=>{t.r(i),t.d(i,{assets:()=>l,contentTitle:()=>s,default:()=>c,frontMatter:()=>r,metadata:()=>o,toc:()=>p});var n=t(5893),a=t(1151);const r={title:"8_PartitionedHeritability",description:"Visualises the partitioned heritability (enrichment)",sidebar_position:3},s="8_PartitionedHeritability",o={id:"ChromOptimise/Pipeline-Explanation/LDSC/PartitionedHeritability",title:"8_PartitionedHeritability",description:"Visualises the partitioned heritability (enrichment)",source:"@site/docs/ChromOptimise/Pipeline-Explanation/LDSC/8_PartitionedHeritability.md",sourceDirName:"ChromOptimise/Pipeline-Explanation/LDSC",slug:"/ChromOptimise/Pipeline-Explanation/LDSC/PartitionedHeritability",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/PartitionedHeritability",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:3,frontMatter:{title:"8_PartitionedHeritability",description:"Visualises the partitioned heritability (enrichment)",sidebar_position:3},sidebar:"documentationSidebar",previous:{title:"SNPAssignment",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/SNPAssignment"},next:{title:"HeritabilityPlots",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots"}},l={},p=[{value:"Explanation",id:"explanation",level:2},{value:"Usage",id:"usage",level:2},{value:"Example usage",id:"example-usage",level:2}];function h(e){const i={a:"a",code:"code",h1:"h1",h2:"h2",p:"p",pre:"pre",...(0,a.a)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(i.h1,{id:"8_partitionedheritability",children:"8_PartitionedHeritability"}),"\n",(0,n.jsx)(i.h2,{id:"explanation",children:"Explanation"}),"\n",(0,n.jsxs)(i.p,{children:["This script will run LDSC again, but now with the ",(0,n.jsx)(i.code,{children:"-h2"})," flag. This signifies\nthat partitioned heritability is wanted. Specifically, what we want to extract\nis the enrichment for each of our states (and the other baseline annotations).\nEnrichment is described as the proportion of heritability divided by the\nproportion of SNPs in the annotated regions in\n",(0,n.jsx)(i.a,{href:"https://www.biorxiv.org/content/10.1101/014241v1.full.pdf",children:"this paper"}),"."]}),"\n",(0,n.jsxs)(i.p,{children:["After extracting this information, the script will call the Rscript\n",(0,n.jsx)(i.a,{href:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots",children:"HeritabilityPlots.R"})," to help with visualisation."]}),"\n",(0,n.jsx)(i.h2,{id:"usage",children:"Usage"}),"\n",(0,n.jsxs)(i.p,{children:["Typically, this script will be ran after being called by\n",(0,n.jsx)(i.a,{href:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/ReferenceLDSCore",children:"7_ReferenceLDSCore.sh"}),". However, one may want to use\nthis script in isolation (perhaps for testing a different set of gwas traits)."]}),"\n",(0,n.jsx)(i.p,{children:"The reason why this script is split away from the previous script in the\npipeline is because the previous script is ran as an array. The memory assigned\nto each task in the array is split evenly and is not dynamic. As a result,\nunless your machine/system has RAM in the 100s of GBs, the partitioned\nheritability step will return out of memory errors. As we only need to compute\nthe plots once, we just use a separate script to run these analyses."}),"\n",(0,n.jsx)(i.h2,{id:"example-usage",children:"Example usage"}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-bash",children:'# Generates heatmap of enrichments and bar plots for all GWAS-traits for the \n# model in the BinSize_200_SampleSize_75_6 directory\nsbatch 7_RunLDSC.sh \\\n--config="path/to/configuration/directory" \\\n--binsize=200 \\\n--samplesize=75 \\\n--nummodels=6 \\\n'})})]})}function c(e={}){const{wrapper:i}={...(0,a.a)(),...e.components};return i?(0,n.jsx)(i,{...e,children:(0,n.jsx)(h,{...e})}):h(e)}},1151:(e,i,t)=>{t.d(i,{Z:()=>o,a:()=>s});var n=t(7294);const a={},r=n.createContext(a);function s(e){const i=n.useContext(r);return n.useMemo((function(){return"function"==typeof e?e(i):{...i,...e}}),[i,e])}function o(e){let i;return i=e.disableParentContext?"function"==typeof e.components?e.components(a):e.components||a:s(e.components),n.createElement(r.Provider,{value:i},e.children)}}}]);