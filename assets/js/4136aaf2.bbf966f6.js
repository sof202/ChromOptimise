"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[9637],{3223:(e,n,i)=>{i.r(n),i.d(n,{assets:()=>r,contentTitle:()=>s,default:()=>c,frontMatter:()=>a,metadata:()=>l,toc:()=>d});var o=i(5893),t=i(1151);const a={title:"0_EGADownloading",description:"The script used to download bam files from EGA",sidebar_position:2},s="0_EGADownloading",l={id:"ChromOptimise/Pipeline-Explanation/EGADownloading",title:"0_EGADownloading",description:"The script used to download bam files from EGA",source:"@site/docs/ChromOptimise/Pipeline-Explanation/0_EGADownloading.md",sourceDirName:"ChromOptimise/Pipeline-Explanation",slug:"/ChromOptimise/Pipeline-Explanation/EGADownloading",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/EGADownloading",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:2,frontMatter:{title:"0_EGADownloading",description:"The script used to download bam files from EGA",sidebar_position:2},sidebar:"documentationSidebar",previous:{title:"Overall explanation of pipeline",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/Overall"},next:{title:"1_MoveFilesToSingleDirectory",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/MoveFilesToSingleDirectory"}},r={},d=[{value:"Explanation",id:"explanation",level:2},{value:"Prerequisites",id:"prerequisites",level:2},{value:"Step one: Login credentials need to be saved in <code>.json</code> format:",id:"step-one-login-credentials-need-to-be-saved-in-json-format",level:3},{value:"Step Two: Create a conda environment with pyega3 installed:",id:"step-two-create-a-conda-environment-with-pyega3-installed",level:3},{value:"Example usage",id:"example-usage",level:2}];function p(e){const n={a:"a",code:"code",h1:"h1",h2:"h2",h3:"h3",p:"p",pre:"pre",...(0,t.a)(),...e.components};return(0,o.jsxs)(o.Fragment,{children:[(0,o.jsx)(n.h1,{id:"0_egadownloading",children:"0_EGADownloading"}),"\n",(0,o.jsx)(n.h2,{id:"explanation",children:"Explanation"}),"\n",(0,o.jsx)(n.p,{children:"This is an artefact of the pipeline originally being built for the blueprint data from EGA."}),"\n",(0,o.jsx)(n.p,{children:"The script will download files from a file of file names using the pyega3 python package. It is likely that this script is not relevant to you, and you can safely ignore it."}),"\n",(0,o.jsx)(n.h2,{id:"prerequisites",children:"Prerequisites"}),"\n",(0,o.jsxs)(n.h3,{id:"step-one-login-credentials-need-to-be-saved-in-json-format",children:["Step one: Login credentials need to be saved in ",(0,o.jsx)(n.code,{children:".json"})," format:"]}),"\n",(0,o.jsx)(n.p,{children:"These should be provided to you via EGA."}),"\n",(0,o.jsx)(n.pre,{children:(0,o.jsx)(n.code,{className:"language-json",children:'{\n    "username":"username@domain.com",\n    "password":"password"\n}\n'})}),"\n",(0,o.jsx)(n.h3,{id:"step-two-create-a-conda-environment-with-pyega3-installed",children:"Step Two: Create a conda environment with pyega3 installed:"}),"\n",(0,o.jsxs)(n.p,{children:["Note down the location of this conda environment. It is to be used in the configuration file ",(0,o.jsx)(n.a,{href:"/ChromOptimise/ChromOptimise/Configuration-Files-Setup#filepathstxt",children:"FilePaths.txt"}),"."]}),"\n",(0,o.jsx)(n.pre,{children:(0,o.jsx)(n.code,{className:"language-bash",children:"# Change 'myenv' to a memorable name\nconda create -n myenv pyega3\n"})}),"\n",(0,o.jsx)(n.h2,{id:"example-usage",children:"Example usage"}),"\n",(0,o.jsx)(n.pre,{children:(0,o.jsx)(n.code,{className:"language-bash",children:'# Downloads all files listed in FileOfFileNames.txt from EGA\nsbatch 0_EGADownloading.sh \\\n--config="path/to/configuration/directory" \\\n--file="FileOfFileNames.txt"\n'})})]})}function c(e={}){const{wrapper:n}={...(0,t.a)(),...e.components};return n?(0,o.jsx)(n,{...e,children:(0,o.jsx)(p,{...e})}):p(e)}},1151:(e,n,i)=>{i.d(n,{Z:()=>l,a:()=>s});var o=i(7294);const t={},a=o.createContext(t);function s(e){const n=o.useContext(a);return o.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function l(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(t):e.components||t:s(e.components),o.createElement(a.Provider,{value:n},e.children)}}}]);