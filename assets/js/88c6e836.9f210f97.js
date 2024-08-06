"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[3888],{7716:(e,s,a)=>{a.r(s),a.d(s,{assets:()=>m,contentTitle:()=>l,default:()=>o,frontMatter:()=>i,metadata:()=>r,toc:()=>h});var n=a(5893),t=a(1151);const i={title:"CalculateBIC",description:"The script used to plot the relative BICs of the models.",sidebar_position:6},l="CalculateBIC",r={id:"ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/CalculateBIC",title:"CalculateBIC",description:"The script used to plot the relative BICs of the models.",source:"@site/docs/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/CalculateBIC.md",sourceDirName:"ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates",slug:"/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/CalculateBIC",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/CalculateBIC",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:6,frontMatter:{title:"CalculateBIC",description:"The script used to plot the relative BICs of the models.",sidebar_position:6},sidebar:"documentationSidebar",previous:{title:"RedundantStateChecker",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/RedundantStateChecker"},next:{title:"PlotLikelihoods",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/PlotLikelihoods"}},m={},h=[{value:"Explanaiton",id:"explanaiton",level:2},{value:"Definition",id:"definition",level:2},{value:"Relative BIC",id:"relative-bic",level:2}];function c(e){const s={annotation:"annotation",h1:"h1",h2:"h2",math:"math",mi:"mi",mn:"mn",mo:"mo",mrow:"mrow",p:"p",semantics:"semantics",span:"span",...(0,t.a)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(s.h1,{id:"calculatebic",children:"CalculateBIC"}),"\n",(0,n.jsx)(s.h2,{id:"explanaiton",children:"Explanaiton"}),"\n",(0,n.jsx)(s.p,{children:'The Bayesian information critereon (BIC) is a heuristic for measuring the\n"information" captured by a model. In contrast with the Akaike information\ncriterion (AIC), BIC is much more punishing on the number of parameters in the\nmodel when the dataset is very large (which it will be in the case of genomic\ndatasets). In general AIC is easier to interpret. However, with large datasets,\nincreasing the number of states is almost always improves the value of this\nheuristic (as the estimated log likelihood is so low).'}),"\n",(0,n.jsx)(s.h2,{id:"definition",children:"Definition"}),"\n",(0,n.jsx)(s.p,{children:"The Bayesian information critereon is given by:"}),"\n",(0,n.jsx)(s.span,{className:"katex-display",children:(0,n.jsxs)(s.span,{className:"katex",children:[(0,n.jsx)(s.span,{className:"katex-mathml",children:(0,n.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",display:"block",children:(0,n.jsxs)(s.semantics,{children:[(0,n.jsxs)(s.mrow,{children:[(0,n.jsx)(s.mi,{children:"B"}),(0,n.jsx)(s.mi,{children:"I"}),(0,n.jsx)(s.mi,{children:"C"}),(0,n.jsx)(s.mo,{children:"="}),(0,n.jsx)(s.mi,{children:"l"}),(0,n.jsx)(s.mi,{children:"n"}),(0,n.jsx)(s.mo,{stretchy:"false",children:"("}),(0,n.jsx)(s.mi,{children:"n"}),(0,n.jsx)(s.mo,{stretchy:"false",children:")"}),(0,n.jsx)(s.mi,{children:"k"}),(0,n.jsx)(s.mo,{children:"\u2212"}),(0,n.jsx)(s.mn,{children:"2"}),(0,n.jsx)(s.mi,{children:"l"}),(0,n.jsx)(s.mi,{children:"n"}),(0,n.jsx)(s.mo,{stretchy:"false",children:"("}),(0,n.jsx)(s.mi,{children:"L"}),(0,n.jsx)(s.mo,{stretchy:"false",children:")"})]}),(0,n.jsx)(s.annotation,{encoding:"application/x-tex",children:"BIC = ln(n)k - 2ln(L)"})]})})}),(0,n.jsxs)(s.span,{className:"katex-html","aria-hidden":"true",children:[(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"0.6833em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.05017em"},children:"B"}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.07847em"},children:"I"}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.07153em"},children:"C"}),(0,n.jsx)(s.span,{className:"mspace",style:{marginRight:"0.2778em"}}),(0,n.jsx)(s.span,{className:"mrel",children:"="}),(0,n.jsx)(s.span,{className:"mspace",style:{marginRight:"0.2778em"}})]}),(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.01968em"},children:"l"}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"n"}),(0,n.jsx)(s.span,{className:"mopen",children:"("}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"n"}),(0,n.jsx)(s.span,{className:"mclose",children:")"}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.03148em"},children:"k"}),(0,n.jsx)(s.span,{className:"mspace",style:{marginRight:"0.2222em"}}),(0,n.jsx)(s.span,{className:"mbin",children:"\u2212"}),(0,n.jsx)(s.span,{className:"mspace",style:{marginRight:"0.2222em"}})]}),(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),(0,n.jsx)(s.span,{className:"mord",children:"2"}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.01968em"},children:"l"}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"n"}),(0,n.jsx)(s.span,{className:"mopen",children:"("}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"L"}),(0,n.jsx)(s.span,{className:"mclose",children:")"})]})]})]})}),"\n",(0,n.jsxs)(s.p,{children:["Where ",(0,n.jsxs)(s.span,{className:"katex",children:[(0,n.jsx)(s.span,{className:"katex-mathml",children:(0,n.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,n.jsxs)(s.semantics,{children:[(0,n.jsx)(s.mrow,{children:(0,n.jsx)(s.mi,{children:"n"})}),(0,n.jsx)(s.annotation,{encoding:"application/x-tex",children:"n"})]})})}),(0,n.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"0.4306em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"n"})]})})]})," is the total number of observations, ",(0,n.jsxs)(s.span,{className:"katex",children:[(0,n.jsx)(s.span,{className:"katex-mathml",children:(0,n.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,n.jsxs)(s.semantics,{children:[(0,n.jsx)(s.mrow,{children:(0,n.jsx)(s.mi,{children:"k"})}),(0,n.jsx)(s.annotation,{encoding:"application/x-tex",children:"k"})]})})}),(0,n.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"0.6944em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",style:{marginRight:"0.03148em"},children:"k"})]})})]})," is the number of parameters\nin the model and ",(0,n.jsxs)(s.span,{className:"katex",children:[(0,n.jsx)(s.span,{className:"katex-mathml",children:(0,n.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,n.jsxs)(s.semantics,{children:[(0,n.jsx)(s.mrow,{children:(0,n.jsx)(s.mi,{children:"L"})}),(0,n.jsx)(s.annotation,{encoding:"application/x-tex",children:"L"})]})})}),(0,n.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"0.6833em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"L"})]})})]})," is the estimated likelihood."]}),"\n",(0,n.jsx)(s.h2,{id:"relative-bic",children:"Relative BIC"}),"\n",(0,n.jsxs)(s.p,{children:["These BIC values are genrally very large and hard to compare. As a result the\nscript will find the BIC of each model relative to the minimum BIC of all\nmodels (this time using a simple ratio as there is no merit to using ",(0,n.jsxs)(s.span,{className:"katex",children:[(0,n.jsx)(s.span,{className:"katex-mathml",children:(0,n.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,n.jsxs)(s.semantics,{children:[(0,n.jsxs)(s.mrow,{children:[(0,n.jsx)(s.mi,{children:"e"}),(0,n.jsx)(s.mi,{children:"x"}),(0,n.jsx)(s.mi,{children:"p"})]}),(0,n.jsx)(s.annotation,{encoding:"application/x-tex",children:"exp"})]})})}),(0,n.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,n.jsxs)(s.span,{className:"base",children:[(0,n.jsx)(s.span,{className:"strut",style:{height:"0.625em",verticalAlign:"-0.1944em"}}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"e"}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"x"}),(0,n.jsx)(s.span,{className:"mord mathnormal",children:"p"})]})})]})," with\nthis heuristic)."]}),"\n",(0,n.jsx)(s.p,{children:"A smaller value for BIC generally indicates a better model. It is important to\nnote that BIC is a heuristics, not a metric. It merely approximates how good a\nmodel is relative to its complexity. As such, the plot and .csv files that come\nfrom this script is only given to provide further information to the user. It\nisn't directly implemented into the pipeline as it isn't concrete."}),"\n",(0,n.jsx)(s.p,{children:"You might use BIC in the following scenario: Perhaps a model with 6 states is\nbetter under this metric than a model with 8 states, despite the more complex\nmodel having no 'redundant states'. Considering the model with 6 states has a\nlower BIC than the model with 8, the user may wish to save on computing power\nand use the 6 state model in further downstream analysis."})]})}function o(e={}){const{wrapper:s}={...(0,t.a)(),...e.components};return s?(0,n.jsx)(s,{...e,children:(0,n.jsx)(c,{...e})}):c(e)}},1151:(e,s,a)=>{a.d(s,{Z:()=>r,a:()=>l});var n=a(7294);const t={},i=n.createContext(t);function l(e){const s=n.useContext(i);return n.useMemo((function(){return"function"==typeof e?e(s):{...s,...e}}),[s,e])}function r(e){let s;return s=e.disableParentContext?"function"==typeof e.components?e.components(t):e.components||t:l(e.components),n.createElement(i.Provider,{value:s},e.children)}}}]);