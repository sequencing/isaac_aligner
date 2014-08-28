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
 ** \file LaneSummary.xsl
 **
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:str="http://exslt.org/strings"
xmlns:math="http://exslt.org/math"
xmlns:isaac="http://www.illumina.com/isaac"
>

<xsl:template name="generateAlignmentStatsCells">
    <xsl:param name="read"/>
    <xsl:variable name="clustersRaw" select="sum(Tile/Raw/ClusterCount)"/>
    <xsl:variable name="clustersPF" select="sum(Tile/Pf/ClusterCount)"/>
    <xsl:variable name="uniquelyAlignedFragmentsPF" select="sum(Tile/Pf/Read[@number=$read]/UniquelyAlignedFragments/Count)"/>
    <xsl:variable name="alignedFragmentsPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/AlignedCount)"/>
        
    <!-- <xsl:variable name="alignScoreSumPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/AlignmentScoreSum)"/> -->

    <xsl:variable name="uniqMismatchesRaw" select="sum(Tile/Raw/Read[@number=$read]/UniquelyAlignedFragments/Mismatches)"/>
    <xsl:variable name="uniqMismatchesPF" select="sum(Tile/Pf/Read[@number=$read]/UniquelyAlignedFragments/Mismatches)"/>
    <xsl:variable name="uniquelyAlignedBasesOutsideIndelsRaw" select="sum(Tile/Raw/Read[@number=$read]/UniquelyAlignedFragments/BasesOutsideIndels)"/>
    <xsl:variable name="uniquelyAlignedBasesOutsideIndelsPF" select="sum(Tile/Pf/Read[@number=$read]/UniquelyAlignedFragments/BasesOutsideIndels)"/>

    <xsl:variable name="allMismatchesRaw" select="sum(Tile/Raw/Read[@number=$read]/AllFragments/Mismatches)"/>
    <xsl:variable name="allMismatchesPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/Mismatches)"/>
    <xsl:variable name="allBasesOutsideIndelsRaw" select="sum(Tile/Raw/Read[@number=$read]/AllFragments/BasesOutsideIndels)"/>
    <xsl:variable name="allBasesOutsideIndelsPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/BasesOutsideIndels)"/>

    <xsl:variable name="yieldPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/Yield)"/>

    <xsl:variable name="yieldPFQ30" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/YieldQ30)"/>

    <xsl:variable name="qualityScoreSumPF" select="sum(Tile/Pf/Read[@number=$read]/AllFragments/QualityScoreSum)"/>
    
    <!--xsl:variable name="referenceName" select="$flowcellNode/Barcode[@number=$barcodeId]/Lane[@number=$lane]/ReferenceName"/-->
    
    <td><xsl:if test="0 != $uniquelyAlignedBasesOutsideIndelsRaw">
        <xsl:value-of select="format-number($uniqMismatchesRaw div $uniquelyAlignedBasesOutsideIndelsRaw * 100, '0.00')"/>/<xsl:value-of select="format-number($allMismatchesRaw div $allBasesOutsideIndelsRaw * 100, '0.00')"/>
    </xsl:if></td>
    <td><xsl:value-of select="format-number($clustersPF, '###,###,###,###,###')"/></td>
    <td><xsl:value-of select="format-number(round($yieldPF div 1000000), '###,###,###,###,###')"/></td>
    <td><xsl:if test="0 != $clustersRaw"><xsl:value-of select="format-number($clustersPF div $clustersRaw * 100, '0.00')"/></xsl:if></td>
    <td><xsl:if test="0 != $clustersPF">
        <xsl:value-of select="format-number($uniquelyAlignedFragmentsPF div $clustersPF * 100, '0.00')"/>/<xsl:value-of select="format-number($alignedFragmentsPF div $clustersPF * 100, '0.00')"/>
    </xsl:if></td>
    <!-- td><xsl:if test="0 != $clustersPF"><xsl:value-of select="format-number($alignScoreSumPF div $clustersPF, '0.00')"/></xsl:if></td-->
    <td><xsl:if test="0 != $uniquelyAlignedBasesOutsideIndelsPF">
        <xsl:value-of select="format-number($uniqMismatchesPF div $uniquelyAlignedBasesOutsideIndelsPF * 100, '0.00')"/>/<xsl:value-of select="format-number($allMismatchesPF div $allBasesOutsideIndelsPF * 100, '0.00')"/>
    </xsl:if></td>
    <td><xsl:if test="0 != $yieldPF"><xsl:value-of select="format-number($yieldPFQ30 div $yieldPF * 100, '0.00')"/></xsl:if></td>
    <td><xsl:if test="0 != $yieldPF"><xsl:value-of select="format-number($qualityScoreSumPF div $yieldPF, '0.00')"/></xsl:if></td>

</xsl:template>

<xsl:template name="generateDemultiplexingStatsCells">
    <xsl:param name="flowcellNode" select="../../../.."/>
    <xsl:param name="projectId" select="../../../@name"/>
    <xsl:param name="sampleId" select="../../@name"/>
    <xsl:param name="barcodeId" select="../@name"/>
    <xsl:variable name="laneNumber" select="@number"/>

    <xsl:variable name="demuxStatsFlowcellNode" 
        select="$DEMULTIPLEXING_STATS_XML/Stats/Flowcell[@flowcell-id=$flowcellNode/@flowcell-id]"/>

    <xsl:variable name="demuxStatsBarcodeLaneNode" select="$demuxStatsFlowcellNode/
Project[@name=$projectId]/Sample[@name=$sampleId]/Barcode[@name=$barcodeId]/Lane[@number=$laneNumber]"/>

    <xsl:variable name="demuxStatsLaneNode" select="$demuxStatsFlowcellNode/
Project[@name='all']/Sample[@name='all']/Barcode[@name='all']/Lane[@number=$laneNumber]"/>

    <xsl:variable name="laneClustersRaw" select="$demuxStatsLaneNode/BarcodeCount"/>
    <xsl:variable name="clustersRaw" select="$demuxStatsBarcodeLaneNode/BarcodeCount"/>
    <xsl:variable name="perfectBarcodeRaw" select="$demuxStatsBarcodeLaneNode/PerfectBarcodeCount"/>
    <xsl:variable name="oneMismatchBarcodeRaw" select="$demuxStatsBarcodeLaneNode/OneMismatchBarcodeCount"/>
    
    <td><xsl:value-of select="format-number($clustersRaw, '###,###,###,###,###')"/></td>
    <td><xsl:if test="0 != $laneClustersRaw"><xsl:value-of select="format-number($clustersRaw div $laneClustersRaw * 100, '0.00')"/></xsl:if></td>
    <td><xsl:if test="0 != $clustersRaw"><xsl:value-of select="format-number($perfectBarcodeRaw div $clustersRaw * 100, '0.00')"/></xsl:if></td>
    <td><xsl:if test="0 != $clustersRaw"><xsl:value-of select="format-number($oneMismatchBarcodeRaw div $clustersRaw * 100, '0.00')"/></xsl:if></td>

</xsl:template>

<xsl:template name="generateExpandedLaneResultsSummaryTable">
    <xsl:param name="flowcellNode" select="../../.."/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    <xsl:param name="read"/>
    <xsl:param name="showBarcodes"/>

    <table border="1" ID="ReportTable">

    <tr>
        <xsl:if test="$showBarcodes and $projectId='all'">
            <th colspan="4">Lane</th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $projectId!='all' and $sampleId='all'">
            <th colspan="3">Lane</th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $projectId!='all' and $sampleId!='all' and $barcodeId='all'">
            <th colspan="2">Lane</th>
        </xsl:if>
        <xsl:if test="true()!=$showBarcodes or ($projectId!='all' and $sampleId!='all' and $barcodeId!='all')">
            <th colspan="1">Lane</th>
        </xsl:if>
        <th colspan="5">Raw data</th>
        <th colspan="7">Filtered data</th>
    </tr>
    <tr>
        <th>#</th>
        <xsl:if test="$showBarcodes and $projectId='all'">
            <th>Project</th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $sampleId='all'">
            <th>Sample</th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $barcodeId='all'">
            <th>Barcode sequence</th>
        </xsl:if>
        <th>Clusters</th>
        <th>% of the<br/>lane</th>
        <th>% Perfect<br/>barcode</th>
        <th>% One mismatch<br/>barcode</th>
        <th>Mismatch % (mapq&gt;3/all)</th>
        <th>Clusters</th>
        <th>Yield (Mbases)</th>
        <th>% PF<br/>Clusters</th>
        <th>% Align (mapq&gt;3/all)</th>
        <!--th>Alignment Score</th-->
        <th>Mismatch % (mapq&gt;3/all)</th>
        <th>% >= Q30<br/>bases</th>
        <th>Mean Quality<br/>Score</th>
    </tr>

    <xsl:for-each select="$flowcellNode/
Project[(true()=$showBarcodes and (($projectId='all' and @name!='all') or ($projectId!='all' and @name=$projectId))) or 
(false()=$showBarcodes and @name=$projectId)]/
Sample[(true()=$showBarcodes and (($sampleId='all' and @name!='all') or ($sampleId!='all' and @name=$sampleId))) or 
(false()=$showBarcodes and @name=$sampleId)]/
Barcode[(true()=$showBarcodes and (($barcodeId='all' and @name!='all') or ($barcodeId!='all' and @name=$barcodeId))) or
(false()=$showBarcodes and @name=$barcodeId)]/
Lane">
    <xsl:sort select="@number"/>

    <tr>
        <td><xsl:value-of select="@number"/></td>
        <xsl:if test="$showBarcodes and $projectId='all'">
            <td><xsl:value-of select="../../../@name"/></td>
        </xsl:if>
        <xsl:if test="$showBarcodes and $sampleId='all'">
            <td><xsl:value-of select="../../@name"/></td>
        </xsl:if>
        <xsl:if test="$showBarcodes and $barcodeId='all'">
            <td><xsl:value-of select="../@name"/></td>
        </xsl:if>
        <xsl:call-template name="generateDemultiplexingStatsCells"/>
        <xsl:call-template name="generateAlignmentStatsCells">
            <xsl:with-param name="read" select="$read"/>
        </xsl:call-template>
    </tr>
    </xsl:for-each>

    </table>

</xsl:template>

<xsl:template name="generateLanePairedStatsTable">
    <xsl:param name="flowcellNode" select="../../.."/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    <xsl:param name="showBarcodes"/>

    <table border="1" ID="ReportTable">
    <xsl:call-template name="generatePairStatsTableHeader">
        <xsl:with-param name="showBarcodes" select="$showBarcodes"/>
        <xsl:with-param name="tableTypeName" select="'Lane'"/>
    </xsl:call-template>

<!--     <xsl:for-each select="$projectNode/Sample[@name=$sampleId]/Barcode[@name=$barcodeId]/Lane">
 -->    <xsl:for-each select="$flowcellNode/
Project[(true()=$showBarcodes and (($projectId='all' and @name!='all') or ($projectId!='all' and @name=$projectId))) or 
(false()=$showBarcodes and @name=$projectId)]/
Sample[(true()=$showBarcodes and (($sampleId='all' and @name!='all') or ($sampleId!='all' and @name=$sampleId))) or 
(false()=$showBarcodes and @name=$sampleId)]/
Barcode[(true()=$showBarcodes and (($barcodeId='all' and @name!='all') or ($barcodeId!='all' and @name=$barcodeId))) or
(false()=$showBarcodes and @name=$barcodeId)]/
Lane">
        <xsl:sort select="@number"/>
        <xsl:variable name="lane" select="@number"/>

        <xsl:call-template name="generatePairStatsTableRow">
            <xsl:with-param name="showProject" select="$showBarcodes and $projectId='all'"/>
            <xsl:with-param name="showSample" select="$showBarcodes and $sampleId='all'"/>
            <xsl:with-param name="showBarcode" select="$showBarcodes and $barcodeId='all'"/>
        </xsl:call-template>
    </xsl:for-each>

    </table>

</xsl:template>


<xsl:template name="generateUnknownBarcodesTable">
    <xsl:param name="flowcellNode" select="../../.."/>

<table border="1" ID="ReportTable">


<tr>
    <th>Lane</th>
    <th>Count</th>
    <th>Sequence</th>
</tr>
<!--
 -->
    <xsl:for-each select="$DEMULTIPLEXING_STATS_XML/Stats/Flowcell[@flowcell-id=$flowcellNode/@flowcell-id]/Lane">
<tr>
    <xsl:element name="th">
        <xsl:attribute name="rowspan"><xsl:value-of select="count(TopUnknownBarcodes/Barcode)"/></xsl:attribute>
        <xsl:value-of select="@number"/>
    </xsl:element>
    <td><xsl:value-of select="TopUnknownBarcodes/Barcode[position()=1]/@count"/></td>
    <td><xsl:value-of select="TopUnknownBarcodes/Barcode[position()=1]/@sequence"/></td>
</tr>
        <xsl:for-each select="TopUnknownBarcodes/Barcode[position()!=1]">
<tr>
    <td><xsl:value-of select="@count"/></td>
    <td><xsl:value-of select="@sequence"/></td>
</tr>
        </xsl:for-each>
    </xsl:for-each>

    </table>

</xsl:template>

</xsl:stylesheet>
