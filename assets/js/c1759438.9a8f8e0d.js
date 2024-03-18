"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[1395],{4691:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>l,contentTitle:()=>r,default:()=>d,frontMatter:()=>o,metadata:()=>a,toc:()=>c});var i=n(5893),s=n(1151);const o={title:"7_RunLDSC",description:"The script used to integrate LDSC.",sidebar_position:1},r="7_RunLDSC",a={id:"ChromOptimise/Pipeline-Explanation/LDSC/RunLDSC",title:"7_RunLDSC",description:"The script used to integrate LDSC.",source:"@site/docs/ChromOptimise/Pipeline-Explanation/LDSC/7_RunLDSC.md",sourceDirName:"ChromOptimise/Pipeline-Explanation/LDSC",slug:"/ChromOptimise/Pipeline-Explanation/LDSC/RunLDSC",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/RunLDSC",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:1,frontMatter:{title:"7_RunLDSC",description:"The script used to integrate LDSC.",sidebar_position:1},sidebar:"documentationSidebar",previous:{title:"LDSC integration",permalink:"/ChromOptimise/category/ldsc-integration"},next:{title:"SNPAssignment",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/SNPAssignment"}},l={},c=[{value:"Explanation",id:"explanation",level:2},{value:"Example usage",id:"example-usage",level:2}];function h(e){const t={a:"a",code:"code",h1:"h1",h2:"h2",li:"li",p:"p",pre:"pre",ul:"ul",...(0,s.a)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(t.h1,{id:"7_runldsc",children:"7_RunLDSC"}),"\n",(0,i.jsx)(t.h2,{id:"explanation",children:"Explanation"}),"\n",(0,i.jsxs)(t.p,{children:["This is a script that will run 2 R scripts and also ",(0,i.jsx)(t.a,{href:"https://github.com/bulik/ldsc",children:"LDSC"}),", which outputs information useful for determining if a ChromHMM model has biologically relevant states (not just statistically relevant)."]}),"\n",(0,i.jsx)(t.p,{children:"Note that you still need to input the bin size, sample size and the number of models learned here. The reason for this is because the file structure is designed such that multiple runs of the same dataset can be analysed concurrently."}),"\n",(0,i.jsxs)(t.p,{children:["The script achieves this using ",(0,i.jsx)(t.a,{href:"https://www.nature.com/articles/ng.3404",children:"partitioned heritability"})," and requires reference files in order to work properly. These files can be downloaded from ",(0,i.jsx)(t.a,{href:"https://zenodo.org/records/10515792",children:"this online repository"}),". The reference files do not have to be from the 1000 genomes project (just any large collection of SNPs will do). The required files are:"]}),"\n",(0,i.jsxs)(t.ul,{children:["\n",(0,i.jsx)(t.li,{children:"1000 genomes (or similar) PLINK files (for your genome build)"}),"\n",(0,i.jsx)(t.li,{children:"1000 genomes (or similar) weight files (built from ldsc, best to get these from the online repository)"}),"\n",(0,i.jsx)(t.li,{children:"Collection of GWAS traits in sumstats format (again, best to get these from the online repository)"}),"\n"]}),"\n",(0,i.jsx)(t.h2,{id:"example-usage",children:"Example usage"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'# Generates heatmap of enrichments and bar plots for all GWAS-traits for the \n# model with 5 states in the model directory\nsbatch 7_RunLDSC.sh \\\n--config="path/to/configuration/directory" \\\n--state=5 \\\n--binsize=200 \\\n--samplesize=75 \\\n--nummodels=6 \\\n'})}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'# Generates heatmap of enrichments and bar plots for all height based traits for  \n# the model with 6 states in the model directory\nsbatch 7_RunLDSC.sh \\\n--config="path/to/configuration/directory" \\\n--state=6 \\\n--gwas="height" \\\n--binsize=200 \\\n--samplesize=75 \\\n--nummodels=6 \\\n'})})]})}function d(e={}){const{wrapper:t}={...(0,s.a)(),...e.components};return t?(0,i.jsx)(t,{...e,children:(0,i.jsx)(h,{...e})}):h(e)}},1151:(e,t,n)=>{n.d(t,{Z:()=>a,a:()=>r});var i=n(7294);const s={},o=i.createContext(s);function r(e){const t=i.useContext(o);return i.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function a(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(s):e.components||s:r(e.components),i.createElement(o.Provider,{value:t},e.children)}}}]);