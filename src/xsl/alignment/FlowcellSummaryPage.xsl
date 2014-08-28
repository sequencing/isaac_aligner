<?xml version="1.0"?>
<!--
/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2014 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** BSD 2-Clause License
 **
 ** You should have received a copy of the BSD 2-Clause License
 ** along with this program. If not, see
 ** <https://github.com/sequencing/licenses/>.
 **
 ** \file FlowcellSummaryPage.xsl
 **
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:isaac="http://www.illumina.com/isaac"
xmlns:str="http://exslt.org/strings"
xmlns:math="http://exslt.org/math"
xmlns:exsl="http://exslt.org/common"
extension-element-prefixes="exsl"
exclude-result-prefixes="str math"
> 


<xsl:template name="generateFlowcellYieldSummaryTable">
    <xsl:param name="projectNode" select="../.."/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    
    <xsl:variable name="clustersPF" select="sum($projectNode/Sample[@name=$sampleId]/Barcode[@name=$barcodeId]/Lane/Tile/Pf/ClusterCount)"/>
    <xsl:variable name="clustersRaw" select="sum($projectNode/Sample[@name=$sampleId]/Barcode[@name=$barcodeId]/Lane/Tile/Raw/ClusterCount)"/>
        
    <xsl:variable name="flowcellYieldPf" select="sum($projectNode/Sample[@name=$sampleId]/Barcode[@name=$barcodeId]/Lane/Tile/Pf/Read/AllFragments/Yield)"/>

    <table border="1" ID="ReportTable">
    <tr><th>Clusters (Raw)</th><th>Clusters(PF)</th><th>Yield (MBases)</th></tr>
    <tr>
        <td><xsl:value-of select="format-number($clustersRaw, '###,###,###,###,###')"/></td>
        <td><xsl:value-of select="format-number($clustersPF, '###,###,###,###,###')"/></td>
        <td><xsl:value-of select="format-number(round($flowcellYieldPf) div 1000000, '###,###,###,###,###')"/></td>
    </tr>
    </table>
</xsl:template>

<xsl:template name="generateFlowcellSummaryTables">
    <xsl:param name="flowcellNode" select="../../.."/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    
    <!--table border="1" ID="ReportTable">
        <tr><td>Flow Cell ID</td><td><xsl:value-of select="$flowcellNode/@flowcell-id"/></td></tr>
        <tr><td>Primary Analysis</td><td><xsl:value-of select="'TODO: Primary Analysis Software version'"/></td></tr>
        <tr><td>Secondary Analysis</td><td><xsl:value-of select="'@iSAAC_VERSION_FULL@'"/></td></tr>
    </table-->

    <xsl:call-template name="generateFlowcellYieldSummaryTable"/>

</xsl:template>

<xsl:template name="generateFlowcellSummaryPage">
    <xsl:param name="flowcellNode" select="../../.."/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    <xsl:param name="showBarcodes"/>
    <xsl:param name="alternativeViewRef"/>
    <xsl:param name="alternativeViewText"/>
<html>
<link rel="stylesheet" href="../../../../{$CSS_FILE_NAME}" type="text/css"/>
<body>

    <table width="100%">
        <tr>
            <td>
    <p><xsl:call-template name="getBarcodeDisplayPath"/></p>
            </td>
            <td>
    <xsl:if test="../../../@flowcell-id!='all' and $barcodeId='all'">
                <p align="right">
    <xsl:element name="a">
        <xsl:attribute name="href"><xsl:value-of select="$alternativeViewRef"/></xsl:attribute>
        <xsl:value-of select="$alternativeViewText"/>
    </xsl:element>
                </p>
    </xsl:if>
            </td>
        </tr>
    </table>
    
    <xsl:call-template name="generateFlowcellSummaryTables"/>
    
<p>Lane Summary : Read 1</p>
    <xsl:call-template name="generateExpandedLaneResultsSummaryTable">
        <xsl:with-param name="read" select="'1'"/>
        <xsl:with-param name="showBarcodes" select="$showBarcodes"/>
    </xsl:call-template>

<p>Lane Summary : Read 2</p>
    <xsl:call-template name="generateExpandedLaneResultsSummaryTable">
        <xsl:with-param name="read" select="'2'"/>
        <xsl:with-param name="showBarcodes" select="$showBarcodes"/>
    </xsl:call-template>
    

<xsl:if test="../../../@flowcell-id!='all'">
<p>
Flowcell Tile Mismatch Graphs
    <xsl:element name="a"><xsl:attribute name="href">../../../../<xsl:call-template name="getFlowcellMismatchGraphsPfLocalPath"/></xsl:attribute>Pf</xsl:element>
 / 
    <xsl:element name="a"><xsl:attribute name="href">../../../../<xsl:call-template name="getFlowcellMismatchGraphsRawLocalPath"/></xsl:attribute>Raw</xsl:element>
</p>
<p>
Flowcell Tile Mismatch Curves
    <xsl:element name="a"><xsl:attribute name="href">../../../../<xsl:call-template name="getFlowcellMismatchCurvesPfLocalPath"/></xsl:attribute>Pf</xsl:element>
 / 
    <xsl:element name="a"><xsl:attribute name="href">../../../../<xsl:call-template name="getFlowcellMismatchCurvesRawLocalPath"/></xsl:attribute>Raw</xsl:element>
</p>
</xsl:if>

<p>Additional Paired Statistics</p>
    <xsl:call-template name="generateLanePairedStatsTable">
        <xsl:with-param name="showBarcodes" select="$showBarcodes"/>
    </xsl:call-template>
     
    <xsl:if test="'unknown'=$barcodeId">
<p>Top Unknown Barcodes</p>
        <xsl:call-template name="generateUnknownBarcodesTable"/>
    </xsl:if>
<p>@iSAAC_VERSION_FULL@</p>

</body>
</html>

</xsl:template>

</xsl:stylesheet>
