import React from 'react';
import ComponentCreator from '@docusaurus/ComponentCreator';

export default [
  {
    path: '/ChromOptimise/markdown-page',
    component: ComponentCreator('/ChromOptimise/markdown-page', 'cc7'),
    exact: true
  },
  {
    path: '/ChromOptimise/',
    component: ComponentCreator('/ChromOptimise/', '6fe'),
    routes: [
      {
        path: '/ChromOptimise/',
        component: ComponentCreator('/ChromOptimise/', 'c48'),
        routes: [
          {
            path: '/ChromOptimise/',
            component: ComponentCreator('/ChromOptimise/', 'fca'),
            routes: [
              {
                path: '/ChromOptimise/category/chromoptimise-documentation',
                component: ComponentCreator('/ChromOptimise/category/chromoptimise-documentation', 'aab'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/ChromHMM-overview',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/ChromHMM-overview', 'b9b'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Configuration-Files-Setup',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Configuration-Files-Setup', '3f3'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Factors-that-affect-the-output',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Factors-that-affect-the-output', '153'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Memory-Profiling',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Memory-Profiling', 'e8a'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Pipeline-Explanation',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Pipeline-Explanation', 'fa1'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Processing-Times',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Processing-Times', 'a58'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/SLURM-Workload-Manager-Information',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/SLURM-Workload-Manager-Information', '265'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/ChromOptimise/Supplementary-pipeline-explanation',
                component: ComponentCreator('/ChromOptimise/ChromOptimise/Supplementary-pipeline-explanation', '73c'),
                exact: true,
                sidebar: "documentationSidebar"
              },
              {
                path: '/ChromOptimise/',
                component: ComponentCreator('/ChromOptimise/', 'b09'),
                exact: true,
                sidebar: "documentationSidebar"
              }
            ]
          }
        ]
      }
    ]
  },
  {
    path: '*',
    component: ComponentCreator('*'),
  },
];
