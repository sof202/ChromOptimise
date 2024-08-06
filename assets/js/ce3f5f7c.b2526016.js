"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[3356],{8034:(e,t,i)=>{i.r(t),i.d(t,{assets:()=>l,contentTitle:()=>a,default:()=>c,frontMatter:()=>o,metadata:()=>r,toc:()=>h});var s=i(5893),n=i(1151);const o={title:"Execute_Redundancy_Metrics",description:"The script that produces the redundancy metrics for a single model file.",sidebar_position:3},a="Execute_Redundancy_Metrics",r={id:"ChromOptimise/Supplementary-Pipeline-Explanation/Execute_Redundancy_Metrics",title:"Execute_Redundancy_Metrics",description:"The script that produces the redundancy metrics for a single model file.",source:"@site/docs/ChromOptimise/Supplementary-Pipeline-Explanation/Execute_Redundancy_Metrics.md",sourceDirName:"ChromOptimise/Supplementary-Pipeline-Explanation",slug:"/ChromOptimise/Supplementary-Pipeline-Explanation/Execute_Redundancy_Metrics",permalink:"/ChromOptimise/ChromOptimise/Supplementary-Pipeline-Explanation/Execute_Redundancy_Metrics",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:3,frontMatter:{title:"Execute_Redundancy_Metrics",description:"The script that produces the redundancy metrics for a single model file.",sidebar_position:3},sidebar:"documentationSidebar",previous:{title:"Generate_Big_Model",permalink:"/ChromOptimise/ChromOptimise/Supplementary-Pipeline-Explanation/Generate_Big_Model"},next:{title:"Factors that affect the output",permalink:"/ChromOptimise/ChromOptimise/Factors-that-affect-the-output"}},l={},h=[{value:"Explanation",id:"explanation",level:2},{value:"Obtaining &#39;good&#39; threshold parameters",id:"obtaining-good-threshold-parameters",level:2},{value:"<code>emissions_threshold</code>",id:"emissions_threshold",level:3},{value:"<code>isolation_threshold</code>",id:"isolation_threshold",level:3},{value:"Example Usage",id:"example-usage",level:2}];function d(e){const t={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",h3:"h3",img:"img",p:"p",pre:"pre",...(0,n.a)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(t.h1,{id:"execute_redundancy_metrics",children:"Execute_Redundancy_Metrics"}),"\n",(0,s.jsx)(t.h2,{id:"explanation",children:"Explanation"}),"\n",(0,s.jsxs)(t.p,{children:["This script runs the same R scripts that 6_OptimalNumberOfStates.sh does. Head\nover to the ",(0,s.jsx)(t.a,{href:"/category/optimal-number-of-states",children:"relevant pages"})," for\nexplanations on the individual scripts. The only difference is that a histogram\nand heatmap is generated for the Euclidean distances calculated (as it can be\nrather difficult to interpret these files when there are lots of states)."]}),"\n",(0,s.jsx)(t.h2,{id:"obtaining-good-threshold-parameters",children:"Obtaining 'good' threshold parameters"}),"\n",(0,s.jsx)(t.p,{children:"Regardless of how these thresholds are chosen, they will always have an element\nof subjectivity to them. A threshold is flawed in this way. But hopefully, the\noutputs of this script (that can be found in subdirectories where the model\nfile inputted resides) can assist in choosing the values for these parameters\nmore intelligently."}),"\n",(0,s.jsx)(t.h3,{id:"emissions_threshold",children:(0,s.jsx)(t.code,{children:"emissions_threshold"})}),"\n",(0,s.jsxs)(t.p,{children:["The ",(0,s.jsx)(t.code,{children:"emissions_threshold"})," is the more important of the two thresholds. The\nhistogram especially should assist in deciding a suitable value for this.\nIdeally your histogram should look something like this:"]}),"\n",(0,s.jsx)(t.p,{children:(0,s.jsx)(t.img,{alt:"Euclidean distances histogram",src:i(1850).Z+"",width:"2100",height:"2100"})}),"\n",(0,s.jsx)(t.p,{children:"From this plot, there is a very obvious gap at around x=0.5. Going below this\nvalue and you see a much larger number of state pairs that have a 'low'\nEuclidean distance between them. The existence of the gap is important as it\nallows one to separate the two groups with a hard line. In the worst case\nscenario the histogram will be completely uniform and it will be impossible to\ntell where the cut off point should be placed. The value of 0.5 is of course,\nstill just a suggestion. One could be more or less stringent with the\nthresholds, just make sure that the value chosen has some meaning behind it."}),"\n",(0,s.jsx)(t.h3,{id:"isolation_threshold",children:(0,s.jsx)(t.code,{children:"isolation_threshold"})}),"\n",(0,s.jsxs)(t.p,{children:["The ",(0,s.jsx)(t.code,{children:"isolation_threshold"})," is less important but still has a large impact on\nyour dataset. This value is unlikely to change with the number of marks in the\ndataset and so it will only need to be changed when a new dataset or chromosome\nis being inspected. The isolation scores of the model analysed will be in the\nsame directory that the model file is in. Opening the text file will reveal\nsomething like:"]}),"\n",(0,s.jsx)(t.pre,{children:(0,s.jsx)(t.code,{className:"language-text",metastring:'title="Isolation_Scores_model-[num].txt"',children:'"state" "isolation_score"\n1 0.0491809430882105\n2 4.51568957599878\n3 5.73507906818568\n4 0.801777139523403\n...\n'})}),"\n",(0,s.jsxs)(t.p,{children:["This file is not converted into a plot as the isolation scores can be very\nlarge for some states (if they were only assigned a few times across the\nchromosome inspected). The general strategy is to see any big jumps in the\nisolation score. Common states will have very low isolation scores, whilst\nrarer states will have much higher isolation scores (or none at all). In this\nlight, the process is very similar to deciding on ",(0,s.jsx)(t.code,{children:"emissions_threshold"}),"."]}),"\n",(0,s.jsx)(t.admonition,{title:"consistency",type:"tip",children:(0,s.jsxs)(t.p,{children:["Make sure to be consistent with your threshold parameters. The threshold\nparameters will likely be ill advised if in your actual analysis you: include\nmore marks, use different bin sizes, inspect different chromosomes, use\ndifferent Phred scores (",(0,s.jsx)(t.em,{children:"etc."}),")."]})}),"\n",(0,s.jsx)(t.h2,{id:"example-usage",children:"Example Usage"}),"\n",(0,s.jsx)(t.pre,{children:(0,s.jsx)(t.code,{className:"language-bash",children:"# Generates redundancy metrics  for the model with 50 states \n# (in big models directory).\nsbatch Execute_Redundancy_Metrics.sh \\\npath/to/config/directory \\\n50\n"})})]})}function c(e={}){const{wrapper:t}={...(0,n.a)(),...e.components};return t?(0,s.jsx)(t,{...e,children:(0,s.jsx)(d,{...e})}):d(e)}},1850:(e,t,i)=>{i.d(t,{Z:()=>s});const s=i.p+"assets/images/Euclidean_distances_histogram_model-80-c2d6b1ac3af0dbe5018ccdb981d0e0b2.png"},1151:(e,t,i)=>{i.d(t,{Z:()=>r,a:()=>a});var s=i(7294);const n={},o=s.createContext(n);function a(e){const t=s.useContext(o);return s.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function r(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(n):e.components||n:a(e.components),s.createElement(o.Provider,{value:t},e.children)}}}]);