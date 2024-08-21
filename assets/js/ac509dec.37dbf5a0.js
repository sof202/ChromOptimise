"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[4240],{3138:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>o,contentTitle:()=>r,default:()=>d,frontMatter:()=>a,metadata:()=>l,toc:()=>h});var s=n(5893),i=n(1151);const a={title:"HeritabilityPlots",description:"The script used to create plots of partitioned heritability.",sidebar_position:4},r="HeritabilityPlots",l={id:"ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots",title:"HeritabilityPlots",description:"The script used to create plots of partitioned heritability.",source:"@site/docs/ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots.md",sourceDirName:"ChromOptimise/Pipeline-Explanation/LDSC",slug:"/ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/HeritabilityPlots",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:4,frontMatter:{title:"HeritabilityPlots",description:"The script used to create plots of partitioned heritability.",sidebar_position:4},sidebar:"documentationSidebar",previous:{title:"5_PartitionedHeritability",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/LDSC/PartitionedHeritability"},next:{title:"Supplementary Pipeline - Usage and Explanation",permalink:"/ChromOptimise/category/supplementary-pipeline---usage-and-explanation"}},o={},h=[{value:"Explanation",id:"explanation",level:2},{value:"Interpretation",id:"interpretation",level:2},{value:"p-value thresholds",id:"p-value-thresholds",level:3},{value:"Warning message",id:"warning-message",level:2}];function c(e){const t={a:"a",admonition:"admonition",annotation:"annotation",blockquote:"blockquote",code:"code",h1:"h1",h2:"h2",h3:"h3",li:"li",math:"math",mi:"mi",mn:"mn",mo:"mo",mrow:"mrow",msup:"msup",ol:"ol",p:"p",semantics:"semantics",span:"span",ul:"ul",...(0,i.a)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(t.h1,{id:"heritabilityplots",children:"HeritabilityPlots"}),"\n",(0,s.jsx)(t.h2,{id:"explanation",children:"Explanation"}),"\n",(0,s.jsxs)(t.p,{children:["This script will take in the ",(0,s.jsx)(t.code,{children:".results"})," files that are outputted from ldsc's\npartitioned heritability and visualise them for you. The two output files are:"]}),"\n",(0,s.jsxs)(t.ol,{children:["\n",(0,s.jsx)(t.li,{children:"A heatmap of the enrichment (Proportion of heritability for state /\nproportion of SNPS for state) for each state against each GWAS trait selected"}),"\n",(0,s.jsx)(t.li,{children:"A barplot for the p-values of the enrichment scores for each state for each\nGWAS trait selected"}),"\n"]}),"\n",(0,s.jsx)(t.p,{children:"This is the main output of runLDSC.sh that enables the user to interpret the\nbiological relevance of their states. If lots of states are showing very\nlittle enrichment for any of the traits that you are considering, this\nsuggests that the model is too specific still (or there is a problem with your\ninput data)."}),"\n",(0,s.jsx)(t.h2,{id:"interpretation",children:"Interpretation"}),"\n",(0,s.jsxs)(t.p,{children:["The colour palette of the heatmap is such that higher enrichment values (",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsxs)(t.mrow,{children:[(0,s.jsx)(t.mo,{children:">"}),(0,s.jsx)(t.mn,{children:"1"})]}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:">1"})]})})}),(0,s.jsxs)(t.span,{className:"katex-html","aria-hidden":"true",children:[(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.5782em",verticalAlign:"-0.0391em"}}),(0,s.jsx)(t.span,{className:"mrel",children:">"}),(0,s.jsx)(t.span,{className:"mspace",style:{marginRight:"0.2778em"}})]}),(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.6444em"}}),(0,s.jsx)(t.span,{className:"mord",children:"1"})]})]})]}),")\nare redish and lower enrichment values (",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsxs)(t.mrow,{children:[(0,s.jsx)(t.mo,{children:"<"}),(0,s.jsx)(t.mn,{children:"1"})]}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:"<1"})]})})}),(0,s.jsxs)(t.span,{className:"katex-html","aria-hidden":"true",children:[(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.5782em",verticalAlign:"-0.0391em"}}),(0,s.jsx)(t.span,{className:"mrel",children:"<"}),(0,s.jsx)(t.span,{className:"mspace",style:{marginRight:"0.2778em"}})]}),(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.6444em"}}),(0,s.jsx)(t.span,{className:"mord",children:"1"})]})]})]}),") are yellowish. Any negative\nenrichment or ridiculously high enrichment (",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsxs)(t.mrow,{children:[(0,s.jsx)(t.mo,{children:">"}),(0,s.jsx)(t.mn,{children:"100"})]}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:">100"})]})})}),(0,s.jsxs)(t.span,{className:"katex-html","aria-hidden":"true",children:[(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.5782em",verticalAlign:"-0.0391em"}}),(0,s.jsx)(t.span,{className:"mrel",children:">"}),(0,s.jsx)(t.span,{className:"mspace",style:{marginRight:"0.2778em"}})]}),(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.6444em"}}),(0,s.jsx)(t.span,{className:"mord",children:"100"})]})]})]}),") is given with a grey box. The\nactual values for the enrichment heatmap are also given in the outputted csv\nfile (Enrichments.csv)"]}),"\n",(0,s.jsx)(t.p,{children:"Enrichment is defined as:"}),"\n",(0,s.jsxs)(t.blockquote,{children:["\n",(0,s.jsx)(t.p,{children:"The proportion of SNP heritability explained divided by the proportion of\nSNPs"}),"\n"]}),"\n",(0,s.jsx)(t.p,{children:"In general you should ignore negative enrichment as\nthis usually stems from negative heritability (which is nonsensicle).For\npositive values however, the following interpretations can be made:"}),"\n",(0,s.jsxs)(t.ul,{children:["\n",(0,s.jsx)(t.li,{children:"An enrichment value greater than 1 (redish) implies that the regions defined\nby the annotation explain a high proportion of snp heritability considering\nthey have a relatively small percentage of all snps."}),"\n",(0,s.jsx)(t.li,{children:"An enrichment value less than 1 (yellowish) implies the opposite. Regions\ndefined by the annotaiton generally explain very little snp heritabilty\nconsidering the total percentage of SNPs they account for."}),"\n",(0,s.jsx)(t.li,{children:"An enrichment value close to 1 is what you would see for a category like\nbase. Base is a category that every SNP falls into (to ensure all SNPs are\naccounted for). Base therefore explains 100% of SNP heritability and contains\n100% of all SNPs, resulting in an enrichment of 1. An enrichment of 1 is not\na significant result of partitioned heritability."}),"\n"]}),"\n",(0,s.jsx)(t.h3,{id:"p-value-thresholds",children:"p-value thresholds"}),"\n",(0,s.jsxs)(t.p,{children:["A baseline value of 0.05 was chosen for a p-value threshold. This was arbitrary\nand the user is free to change this in the associated script\n(",(0,s.jsx)(t.code,{children:"HeritabilityPlots.R"}),"). However, this p-value threshold is under jeopardy due\nto the high number of hypotheses that are being tested at once. If you are\nlooking at ~10 different gwas traits in your analysis, you end up testing\nhundreds of hypotheses at once (as you are looking at the enrichment for\nevery annotation as well)."]}),"\n",(0,s.jsx)(t.p,{children:"To account for this, two correction methods were implemented: Bonferroni\ncorrection and FDR correction (using Benjamini and Hochberg method). The FDR\nthreshold is given by the grey dotted line on the bar plots and Bonferroni\nthreshold is given by the black dotted line. These correction methods help\nmitigate the possibility of accepting a low pvalue that ocurred due to chance."}),"\n",(0,s.jsx)(t.p,{children:"The heatmap gives you a quick overview of the p-values with stars. A single\nstar means that the enrichment is fdr significant, 2 stars means the enrichment\nis Bonferroni significant. For more indepth information look at the bar charts\nfor each gwas trait."}),"\n",(0,s.jsxs)(t.admonition,{title:"Number of hypotheses",type:"info",children:[(0,s.jsx)(t.p,{children:"The number of hypotheses should be equal to the number of entries in the\ncomplete (all categories) heatmap. However, the user is not actaully able to\ncontrol the number of annotations that come from the baseline files provided\nby the team behind ldsc. As such, we decided to discount the hypotheses\ngenerated by baseline annotations when calculating corrected p-value\nthresholds. This makes the thresholds more reactive to the user's inputs and\ndata set."}),(0,s.jsx)(t.p,{children:"Also note that these thresholds are to be taken with a mound of salt. Not only\nis the original threshold arbitrary, but the corrected thresholds implicitly\nassume that the hypotheses are independent of one another. In general, gwas\nresults are not independent of one another, so values are over corrected. Keep\nthis in consideration when rejecting the null hypothesis."})]}),"\n",(0,s.jsx)(t.h2,{id:"warning-message",children:"Warning message"}),"\n",(0,s.jsxs)(t.p,{children:["Sometimes, this script will output the file WARNING.txt into the root of the\nplots folder. This file will tell the user that their enrichment heatmap has\na high proportion of negative enrichments. Mathematically a negative enrichment\nis nonsensicle (it is defined as the quotient of two scrictly positive values).\nHowever, ldsc uses an unbiased estimator of ",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsx)(t.mrow,{children:(0,s.jsxs)(t.msup,{children:[(0,s.jsx)(t.mi,{children:"r"}),(0,s.jsx)(t.mn,{children:"2"})]})}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:"r^2"})]})})}),(0,s.jsx)(t.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.8141em"}}),(0,s.jsxs)(t.span,{className:"mord",children:[(0,s.jsx)(t.span,{className:"mord mathnormal",style:{marginRight:"0.02778em"},children:"r"}),(0,s.jsx)(t.span,{className:"msupsub",children:(0,s.jsx)(t.span,{className:"vlist-t",children:(0,s.jsx)(t.span,{className:"vlist-r",children:(0,s.jsx)(t.span,{className:"vlist",style:{height:"0.8141em"},children:(0,s.jsxs)(t.span,{style:{top:"-3.063em",marginRight:"0.05em"},children:[(0,s.jsx)(t.span,{className:"pstrut",style:{height:"2.7em"}}),(0,s.jsx)(t.span,{className:"sizing reset-size6 size3 mtight",children:(0,s.jsx)(t.span,{className:"mord mtight",children:"2"})})]})})})})})]})]})})]})," due to the biases when\nrunning over a large number of SNPS (which can be negative). The creators of\nldsc state that:"]}),"\n",(0,s.jsxs)(t.blockquote,{children:["\n",(0,s.jsx)(t.p,{children:"If more than a few percent of all SNPs have negative LD Scores, this\nprobably means that something is wrong, either that the sample size used for\nestimating LD Scores is too small or the window size used is too large."}),"\n"]}),"\n",(0,s.jsx)(t.p,{children:"If this does happen to you, please look at the enrichment heatmap and check the\nfollowing:"}),"\n",(0,s.jsxs)(t.ul,{children:["\n",(0,s.jsxs)(t.li,{children:["Is there a specific category that appears to be more often negatively\nenriched than not?","\n",(0,s.jsxs)(t.ul,{children:["\n",(0,s.jsxs)(t.li,{children:["If so, then add the category to the list of blacklisted categories in\n",(0,s.jsx)(t.code,{children:"CreateAnnotationFile.R"})," (the vector can be found at the end of the file)."]}),"\n"]}),"\n"]}),"\n",(0,s.jsxs)(t.li,{children:["Are the negative enrichment values somewhat evenly distributed throughout the\nheatmap?","\n",(0,s.jsxs)(t.ul,{children:["\n",(0,s.jsxs)(t.li,{children:["If so, then make the window size smaller in ",(0,s.jsx)(t.code,{children:"7_ReferenceLDSCore.sh"}),". Change\nthe flag ",(0,s.jsx)(t.code,{children:"--ld-wind-cm 1"})," to have a smaller value than 1."]}),"\n",(0,s.jsx)(t.li,{children:"This line can be found in the reference LD scores section of the script"}),"\n"]}),"\n"]}),"\n"]}),"\n",(0,s.jsxs)(t.p,{children:["If problems continue to persist, there may be a problem with the dataset that\nyou are calculating LD scores for. Create an issue for this repository or\n",(0,s.jsx)(t.a,{href:"https://github.com/bulik/ldsc/issues",children:"ldsc"})," for further assistance."]})]})}function d(e={}){const{wrapper:t}={...(0,i.a)(),...e.components};return t?(0,s.jsx)(t,{...e,children:(0,s.jsx)(c,{...e})}):c(e)}},1151:(e,t,n)=>{n.d(t,{Z:()=>l,a:()=>r});var s=n(7294);const i={},a=s.createContext(i);function r(e){const t=s.useContext(a);return s.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function l(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(i):e.components||i:r(e.components),s.createElement(a.Provider,{value:t},e.children)}}}]);