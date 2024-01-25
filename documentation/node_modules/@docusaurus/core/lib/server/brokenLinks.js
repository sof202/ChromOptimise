"use strict";
/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */
Object.defineProperty(exports, "__esModule", { value: true });
exports.handleBrokenLinks = void 0;
const tslib_1 = require("tslib");
const lodash_1 = tslib_1.__importDefault(require("lodash"));
const logger_1 = tslib_1.__importDefault(require("@docusaurus/logger"));
const react_router_config_1 = require("react-router-config");
const utils_1 = require("@docusaurus/utils");
const utils_2 = require("./utils");
function getBrokenLinksForPage({ collectedLinks, pagePath, pageLinks, routes, }) {
    // console.log('routes:', routes);
    function isPathBrokenLink(linkPath) {
        const matchedRoutes = [linkPath.pathname, decodeURI(linkPath.pathname)]
            // @ts-expect-error: React router types RouteConfig with an actual React
            // component, but we load route components with string paths.
            // We don't actually access component here, so it's fine.
            .map((l) => (0, react_router_config_1.matchRoutes)(routes, l))
            .flat();
        return matchedRoutes.length === 0;
    }
    function isAnchorBrokenLink(linkPath) {
        const { pathname, hash } = linkPath;
        // Link has no hash: it can't be a broken anchor link
        if (hash === undefined) {
            return false;
        }
        const targetPage = collectedLinks[pathname] || collectedLinks[decodeURI(pathname)];
        // link with anchor to a page that does not exist (or did not collect any
        // link/anchor) is considered as a broken anchor
        if (!targetPage) {
            return true;
        }
        // it's a broken anchor if the target page exists
        // but the anchor does not exist on that page
        return !targetPage.anchors.includes(hash);
    }
    const brokenLinks = pageLinks.flatMap((link) => {
        const linkPath = (0, utils_1.parseURLPath)(link, pagePath);
        if (isPathBrokenLink(linkPath)) {
            return [
                {
                    link,
                    resolvedLink: (0, utils_1.serializeURLPath)(linkPath),
                    anchor: false,
                },
            ];
        }
        if (isAnchorBrokenLink(linkPath)) {
            return [
                {
                    link,
                    resolvedLink: (0, utils_1.serializeURLPath)(linkPath),
                    anchor: true,
                },
            ];
        }
        return [];
    });
    return brokenLinks;
}
/**
 * The route defs can be recursive, and have a parent match-all route. We don't
 * want to match broken links like /docs/brokenLink against /docs/*. For this
 * reason, we only consider the "final routes" that do not have subroutes.
 * We also need to remove the match-all 404 route
 */
function filterIntermediateRoutes(routesInput) {
    const routesWithout404 = routesInput.filter((route) => route.path !== '*');
    return (0, utils_2.getAllFinalRoutes)(routesWithout404);
}
function getBrokenLinks({ collectedLinks, routes, }) {
    const filteredRoutes = filterIntermediateRoutes(routes);
    return lodash_1.default.mapValues(collectedLinks, (pageCollectedData, pagePath) => getBrokenLinksForPage({
        collectedLinks,
        pageLinks: pageCollectedData.links,
        pageAnchors: pageCollectedData.anchors,
        pagePath,
        routes: filteredRoutes,
    }));
}
function brokenLinkMessage(brokenLink) {
    const showResolvedLink = brokenLink.link !== brokenLink.resolvedLink;
    return `${brokenLink.link}${showResolvedLink ? ` (resolved as: ${brokenLink.resolvedLink})` : ''}`;
}
function createBrokenLinksMessage(pagePath, brokenLinks) {
    const type = brokenLinks[0]?.anchor === true ? 'anchor' : 'link';
    const anchorMessage = brokenLinks.length > 0
        ? `- Broken ${type} on source page path = ${pagePath}:
   -> linking to ${brokenLinks
            .map(brokenLinkMessage)
            .join('\n   -> linking to ')}`
        : '';
    return `${anchorMessage}`;
}
function createBrokenAnchorsMessage(brokenAnchors) {
    if (Object.keys(brokenAnchors).length === 0) {
        return undefined;
    }
    return `Docusaurus found broken anchors!

Please check the pages of your site in the list below, and make sure you don't reference any anchor that does not exist.
Note: it's possible to ignore broken anchors with the 'onBrokenAnchors' Docusaurus configuration, and let the build pass.

Exhaustive list of all broken anchors found:
${Object.entries(brokenAnchors)
        .map(([pagePath, brokenLinks]) => createBrokenLinksMessage(pagePath, brokenLinks))
        .join('\n')}
`;
}
function createBrokenPathsMessage(brokenPathsMap) {
    if (Object.keys(brokenPathsMap).length === 0) {
        return undefined;
    }
    /**
     * If there's a broken link appearing very often, it is probably a broken link
     * on the layout. Add an additional message in such case to help user figure
     * this out. See https://github.com/facebook/docusaurus/issues/3567#issuecomment-706973805
     */
    function getLayoutBrokenLinksHelpMessage() {
        const flatList = Object.entries(brokenPathsMap).flatMap(([pagePage, brokenLinks]) => brokenLinks.map((brokenLink) => ({ pagePage, brokenLink })));
        const countedBrokenLinks = lodash_1.default.countBy(flatList, (item) => item.brokenLink.link);
        const FrequencyThreshold = 5; // Is this a good value?
        const frequentLinks = Object.entries(countedBrokenLinks)
            .filter(([, count]) => count >= FrequencyThreshold)
            .map(([link]) => link);
        if (frequentLinks.length === 0) {
            return '';
        }
        return logger_1.default.interpolate `

It looks like some of the broken links we found appear in many pages of your site.
Maybe those broken links appear on all pages through your site layout?
We recommend that you check your theme configuration for such links (particularly, theme navbar and footer).
Frequent broken links are linking to:${frequentLinks}`;
    }
    return `Docusaurus found broken links!

Please check the pages of your site in the list below, and make sure you don't reference any path that does not exist.
Note: it's possible to ignore broken links with the 'onBrokenLinks' Docusaurus configuration, and let the build pass.${getLayoutBrokenLinksHelpMessage()}

Exhaustive list of all broken links found:
${Object.entries(brokenPathsMap)
        .map(([pagePath, brokenPaths]) => createBrokenLinksMessage(pagePath, brokenPaths))
        .join('\n')}
`;
}
function splitBrokenLinks(brokenLinks) {
    const brokenPaths = {};
    const brokenAnchors = {};
    Object.entries(brokenLinks).forEach(([pathname, pageBrokenLinks]) => {
        const [anchorBrokenLinks, pathBrokenLinks] = lodash_1.default.partition(pageBrokenLinks, (link) => link.anchor);
        if (pathBrokenLinks.length > 0) {
            brokenPaths[pathname] = pathBrokenLinks;
        }
        if (anchorBrokenLinks.length > 0) {
            brokenAnchors[pathname] = anchorBrokenLinks;
        }
    });
    return { brokenPaths, brokenAnchors };
}
function reportBrokenLinks({ brokenLinks, onBrokenLinks, onBrokenAnchors, }) {
    // We need to split the broken links reporting in 2 for better granularity
    // This is because we need to report broken path/anchors independently
    // For v3.x retro-compatibility, we can't throw by default for broken anchors
    // TODO Docusaurus v4: make onBrokenAnchors throw by default?
    const { brokenPaths, brokenAnchors } = splitBrokenLinks(brokenLinks);
    const pathErrorMessage = createBrokenPathsMessage(brokenPaths);
    if (pathErrorMessage) {
        logger_1.default.report(onBrokenLinks)(pathErrorMessage);
    }
    const anchorErrorMessage = createBrokenAnchorsMessage(brokenAnchors);
    if (anchorErrorMessage) {
        logger_1.default.report(onBrokenAnchors)(anchorErrorMessage);
    }
}
async function handleBrokenLinks({ collectedLinks, onBrokenLinks, onBrokenAnchors, routes, }) {
    if (onBrokenLinks === 'ignore' && onBrokenAnchors === 'ignore') {
        return;
    }
    const brokenLinks = getBrokenLinks({ routes, collectedLinks });
    reportBrokenLinks({ brokenLinks, onBrokenLinks, onBrokenAnchors });
}
exports.handleBrokenLinks = handleBrokenLinks;
