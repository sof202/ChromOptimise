"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[7844],{186:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>o,contentTitle:()=>l,default:()=>m,frontMatter:()=>i,metadata:()=>r,toc:()=>h});var s=n(5893),a=n(1151);const i={title:"FlankingStates",description:"The script used to calculate the most likely flanks.",sidebar_position:3},l="FlankingStates",r={id:"ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/FlankingStates",title:"FlankingStates",description:"The script used to calculate the most likely flanks.",source:"@site/docs/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/FlankingStates.md",sourceDirName:"ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates",slug:"/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/FlankingStates",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/FlankingStates",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:3,frontMatter:{title:"FlankingStates",description:"The script used to calculate the most likely flanks.",sidebar_position:3},sidebar:"documentationSidebar",previous:{title:"SimilarEmissions",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/SimilarEmissions"},next:{title:"IsolationScores",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation/OptimalNumberOfStates/IsolationScores"}},o={},h=[{value:"Explanation",id:"explanation",level:2}];function c(e){const t={annotation:"annotation",blockquote:"blockquote",h1:"h1",h2:"h2",math:"math",mi:"mi",mrow:"mrow",p:"p",semantics:"semantics",span:"span",table:"table",tbody:"tbody",td:"td",th:"th",thead:"thead",tr:"tr",...(0,a.a)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(t.h1,{id:"flankingstates",children:"FlankingStates"}),"\n",(0,s.jsx)(t.h2,{id:"explanation",children:"Explanation"}),"\n",(0,s.jsx)(t.p,{children:"The definition of 'flanking state' is as follows:"}),"\n",(0,s.jsxs)(t.blockquote,{children:["\n",(0,s.jsxs)(t.p,{children:["For a bin that is assigned state ",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsx)(t.mrow,{children:(0,s.jsx)(t.mi,{children:"x"})}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:"x"})]})})}),(0,s.jsx)(t.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.4306em"}}),(0,s.jsx)(t.span,{className:"mord mathnormal",children:"x"})]})})]}),", the downstream flanking state is the\nstate assignment for the next bin that is not also assigned state ",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsx)(t.mrow,{children:(0,s.jsx)(t.mi,{children:"x"})}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:"x"})]})})}),(0,s.jsx)(t.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.4306em"}}),(0,s.jsx)(t.span,{className:"mord mathnormal",children:"x"})]})})]}),". The\nupstream flanking state is the same but for the previous bin that is not also\nassigned ",(0,s.jsxs)(t.span,{className:"katex",children:[(0,s.jsx)(t.span,{className:"katex-mathml",children:(0,s.jsx)(t.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(t.semantics,{children:[(0,s.jsx)(t.mrow,{children:(0,s.jsx)(t.mi,{children:"x"})}),(0,s.jsx)(t.annotation,{encoding:"application/x-tex",children:"x"})]})})}),(0,s.jsx)(t.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(t.span,{className:"base",children:[(0,s.jsx)(t.span,{className:"strut",style:{height:"0.4306em"}}),(0,s.jsx)(t.span,{className:"mord mathnormal",children:"x"})]})})]}),"."]}),"\n"]}),"\n",(0,s.jsx)(t.p,{children:"Consider the below example:"}),"\n",(0,s.jsxs)(t.table,{children:[(0,s.jsx)(t.thead,{children:(0,s.jsxs)(t.tr,{children:[(0,s.jsx)(t.th,{children:"Bin number"}),(0,s.jsx)(t.th,{children:"1"}),(0,s.jsx)(t.th,{children:"2"}),(0,s.jsx)(t.th,{children:"3"}),(0,s.jsx)(t.th,{children:"4"}),(0,s.jsx)(t.th,{children:"5"}),(0,s.jsx)(t.th,{children:"6"}),(0,s.jsx)(t.th,{children:"7"})]})}),(0,s.jsx)(t.tbody,{children:(0,s.jsxs)(t.tr,{children:[(0,s.jsx)(t.td,{children:"state assignment"}),(0,s.jsx)(t.td,{children:"1"}),(0,s.jsx)(t.td,{children:"1"}),(0,s.jsx)(t.td,{children:"1"}),(0,s.jsx)(t.td,{children:"2"}),(0,s.jsx)(t.td,{children:"2"}),(0,s.jsx)(t.td,{children:"3"}),(0,s.jsx)(t.td,{children:"2"})]})})]}),"\n",(0,s.jsx)(t.p,{children:"In this example, lets find the flanking states for bin number 4. The next bin\n(bin 5) is also assigned state 2, so we go to the next bin (and so on) until we\nreach a bin that is not assigned 2. In this case, we go to bin 6 and find the\ndownstream flank is 3. In the same way we can go backwards to find the upstream\nflank is 1 (assignment for bin 3)."}),"\n",(0,s.jsx)(t.p,{children:"The most likely flanking states are the expected up/downstream flanking states\nfor any given instance of a state in the state assignment file. These most\nlikely flanking states are what this script aims to find."}),"\n",(0,s.jsx)(t.p,{children:"The most likely flanking states are found using the transition matrix for the\nmodel being inspected. ChromHMM doesn't just factor in the transition\nprobabilities when finding the most likely state assignment. However, using the\nprobabilities from these transition matrices is more generalisable to any set\nof observations (and also saves massively on processing time)."})]})}function m(e={}){const{wrapper:t}={...(0,a.a)(),...e.components};return t?(0,s.jsx)(t,{...e,children:(0,s.jsx)(c,{...e})}):c(e)}},1151:(e,t,n)=>{n.d(t,{Z:()=>r,a:()=>l});var s=n(7294);const a={},i=s.createContext(a);function l(e){const t=s.useContext(i);return s.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function r(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(a):e.components||a:l(e.components),s.createElement(i.Provider,{value:t},e.children)}}}]);