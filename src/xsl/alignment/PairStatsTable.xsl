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
 ** \file PairStatsTable.xsl
 **
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:str="http://exslt.org/strings"
xmlns:math="http://exslt.org/math"
xmlns:casava="http://www.illumina.com/casava/alignment"
> 

<xsl:template name="generatePairStatsTableHeader">
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    <xsl:param name="showBarcodes"/>
    <xsl:param name="tableTypeName"/>
    <tr>
        <xsl:if test="$showBarcodes and $projectId='all'">
            <th colspan="4"><xsl:value-of select="$tableTypeName"/></th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $projectId!='all' and $sampleId='all'">
            <th colspan="3"><xsl:value-of select="$tableTypeName"/></th>
        </xsl:if>
        <xsl:if test="$showBarcodes and $projectId!='all' and $sampleId!='all' and $barcodeId='all'">
            <th colspan="2"><xsl:value-of select="$tableTypeName"/></th>
        </xsl:if>
        <xsl:if test="true()!=$showBarcodes or ($projectId!='all' and $sampleId!='all' and $barcodeId!='all')">
            <th colspan="1"><xsl:value-of select="$tableTypeName"/></th>
        </xsl:if>
        <th rowspan="2">Ref</th>
        <th colspan="5">Relative Orientation Statistics</th>
        <th colspan="5">Mean Template Length Statistics across tiles</th>
        <th colspan="3">Template Statistics<br/>(% of individually uniquely alignable pairs)</th>
    </tr>
    <tr>
        <!--th><xsl:value-of select="$tableTypeName"/></th-->
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
        <th>F-:<br/>&gt;R2 R1&gt;</th>
        <th>F+:<br/>&gt;R1 R2&gt;</th>
        <th>R-:<br/>&lt;R2 R1&gt;</th>
        <th>R+:<br/>&gt;R1 R2&lt;</th>
        <th>Total</th>
        <th>Median</th>
        <th>Below<br/>median SD</th>
        <th>Above<br/>median SD</th>
        <th>Low<br/>thresh.</th>
        <th>High<br/>thresh.</th>
        <th>Too<br/>small</th>
        <th>Too<br/>large</th>
        <th>Orientation<br/>and size OK</th>
    </tr>
</xsl:template>

<xsl:template name="generatePairStatsTableRow">
    <xsl:param name="flowcellNode" select="../../../.."/>
    <xsl:param name="projectId" select="../../../@name"/>
    <xsl:param name="sampleId" select="../../@name"/>
    <xsl:param name="barcodeId" select="../@name"/>
    <xsl:param name="showProject"/>
    <xsl:param name="showSample"/>
    <xsl:param name="showBarcode"/>
    <xsl:param name="lane" select="@number"/>
    
    <xsl:variable name="referenceName" select="$flowcellNode/Barcode[@name=$barcodeId]/Lane[@number=$lane]/ReferenceName"/>
    
    <xsl:variable name="orientationNode" select="Tile/Pf/AlignmentModel"/>
    <xsl:variable name="fm" select="sum($orientationNode/FFp) + sum($orientationNode/RRm)"/>
    <xsl:variable name="fp" select="sum($orientationNode/RRp) + sum($orientationNode/FFm)"/>
    <xsl:variable name="rm" select="sum($orientationNode/RFp) + sum($orientationNode/FRm)"/>
    <xsl:variable name="rp" select="sum($orientationNode/FRp) + sum($orientationNode/RFm)"/>
    <xsl:variable name="orientationTotal" select="$fm + $fp + $rm + $rp"/>

    <xsl:variable name="meanTemplateMedian">
        <xsl:call-template name="meanAndStdev">
            <xsl:with-param name="nodes" select="Tile/AssumedTemplateLength/Median"/>
            <xsl:with-param name="meanFmt" select="'0'"/>
            <xsl:with-param name="devFmt" select="'0'"/>
        </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="meanTemplateLowStdDev">
        <xsl:call-template name="meanAndStdev">
            <xsl:with-param name="nodes" select="Tile/AssumedTemplateLength/LowStdDev"/>
            <xsl:with-param name="meanFmt" select="'0'"/>
            <xsl:with-param name="devFmt" select="'0'"/>
        </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="meanTemplateHighStdDev">
        <xsl:call-template name="meanAndStdev">
            <xsl:with-param name="nodes" select="Tile/AssumedTemplateLength/HighStdDev"/>
            <xsl:with-param name="meanFmt" select="'0'"/>
            <xsl:with-param name="devFmt" select="'0'"/>
        </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="meanTemplateMin">
        <xsl:call-template name="meanAndStdev">
            <xsl:with-param name="nodes" select="Tile/AssumedTemplateLength/Min"/>
            <xsl:with-param name="meanFmt" select="'0'"/>
            <xsl:with-param name="devFmt" select="'0'"/>
        </xsl:call-template>
    </xsl:variable>

    <xsl:variable name="meanTemplateMax">
        <xsl:call-template name="meanAndStdev">
            <xsl:with-param name="nodes" select="Tile/AssumedTemplateLength/Max"/>
            <xsl:with-param name="meanFmt" select="'0'"/>
            <xsl:with-param name="devFmt" select="'0'"/>
        </xsl:call-template>
    </xsl:variable>

    <tr>
        <td><xsl:value-of select="$lane"/></td>
        <xsl:if test="$showProject">
            <td><xsl:value-of select="../../../@name"/></td>
        </xsl:if>
        <xsl:if test="$showSample">
            <td><xsl:value-of select="../../@name"/></td>
        </xsl:if>
        <xsl:if test="$showBarcode">
            <td><xsl:value-of select="../@name"/></td>
        </xsl:if>
        <td><xsl:value-of select="$referenceName"/></td>
        <xsl:if test="0 = $orientationTotal"><td/><td/><td/><td/></xsl:if>
        <xsl:if test="0 != $orientationTotal">
            <td>
                <xsl:value-of select="format-number($fm, '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * $fm div $orientationTotal, '0.00')"/>%)
            </td>
            <td>
                <xsl:value-of select="format-number($fp, '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * $fp div $orientationTotal, '0.00')"/>%)
            </td>
            <td>
                <xsl:value-of select="format-number($rm, '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * $rm div $orientationTotal, '0.00')"/>%)
            </td>
            <td>
                <xsl:value-of select="format-number($rp, '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * $rp div $orientationTotal, '0.00')"/>%)
            </td>
        </xsl:if>
        <td>
            <xsl:value-of select="format-number($orientationTotal, '###,###,###,###,###')"/>
        </td>
        <td><xsl:value-of select="$meanTemplateMedian"/></td>
        <td><xsl:value-of select="$meanTemplateLowStdDev"/></td>
        <td><xsl:value-of select="$meanTemplateHighStdDev"/></td>
        <td><xsl:value-of select="$meanTemplateMin"/></td>
        <td><xsl:value-of select="$meanTemplateMax"/></td>
        <xsl:if test="0 = $orientationTotal"><td/><td/><td/></xsl:if>
        <xsl:if test="0 != $orientationTotal">
            <td>
                <xsl:value-of select="format-number(sum($orientationNode/Undersized), '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * sum($orientationNode/Undersized) div $orientationTotal, '0.00')"/>%)
            </td>
            <td>
                <xsl:value-of select="format-number(sum($orientationNode/Oversized), '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * sum($orientationNode/Oversized) div $orientationTotal, '0.00')"/>%)
            </td>
            <td>
                <xsl:value-of select="format-number(sum($orientationNode/Nominal), '###,###,###,###,###')"/><br/>
                (<xsl:value-of select="format-number(100.0 * sum($orientationNode/Nominal) div $orientationTotal, '0.00')"/>%)
            </td>
        </xsl:if>
    </tr>

</xsl:template>


</xsl:stylesheet>
