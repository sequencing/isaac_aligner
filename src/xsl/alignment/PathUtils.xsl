<?xml version="1.0"?>
<!--
/**
 ** Isaac Genome Alignment Software
 ** Copyright (c) 2010-2012 Illumina, Inc.
 **
 ** This software is provided under the terms and conditions of the
 ** Illumina Open Source Software License 1.
 **
 ** You should have received a copy of the Illumina Open Source
 ** Software License 1 along with this program. If not, see
 ** <https://github.com/downloads/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file PathUtils.xsl
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

<xsl:variable name="MISMATCHES_FILE_SUFFIX" select="'_mismatches.gif'"/>
<xsl:variable name="MISMATCHES_THUMBNAIL_FILE_SUFFIX" select="'_mismatches_thumb.gif'"/>
<xsl:variable name="MISMATCH_CURVES_FILE_SUFFIX" select="'_mismatch-fragments.gif'"/>
<xsl:variable name="MISMATCH_CURVES_THUMBNAIL_FILE_SUFFIX" select="'_mismatch-fragments_thumb.gif'"/>


<xsl:template name="getFlowcellSampleBarcodeSimpleLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/lane.html')"/>
</xsl:template>

<xsl:template name="getFlowcellSampleBarcodeSimpleGlobalPath">
    <xsl:variable name="flowcellSummaryPageFileName" ><xsl:call-template name="getFlowcellSampleBarcodeSimpleLocalPath"/></xsl:variable>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/', $flowcellSummaryPageFileName)"/>
</xsl:template>

<xsl:template name="getFlowcellSampleBarcodeLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="../../@name"/>
    <xsl:param name="sampleId" select="../@name"/>
    <xsl:param name="barcodeId" select="@name"/>
    
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/laneBarcode.html')"/>
</xsl:template>

<xsl:template name="getFlowcellSampleBarcodeGlobalPath">
    <xsl:variable name="flowcellSummaryPageFileName" ><xsl:call-template name="getFlowcellSampleBarcodeLocalPath"/></xsl:variable>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/', $flowcellSummaryPageFileName)"/>
</xsl:template>


<xsl:template name="getFlowcellRefsLocalPath">
    <xsl:value-of select="'tree.html'"/>
</xsl:template>

<xsl:template name="getFlowcellRefsGlobalPath">
    <xsl:variable name="flowcellRefsLocalPath" ><xsl:call-template name="getFlowcellRefsLocalPath"/></xsl:variable>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/', $flowcellRefsLocalPath)"/>
</xsl:template>


<xsl:template name="getFlowcellMismatchGraphsPfLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="'all'"/>
    <xsl:param name="sampleId" select="'all'"/>
    <xsl:param name="barcodeId" select="'all'"/>
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/PfCycleMismatches.html')"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchGraphsRawLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="'all'"/>
    <xsl:param name="sampleId" select="'all'"/>
    <xsl:param name="barcodeId" select="'all'"/>
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/RawCycleMismatches.html')"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchGraphsPfGlobalPath">
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/')"/><xsl:call-template name="getFlowcellMismatchGraphsPfLocalPath"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchGraphsRawGlobalPath">
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/')"/><xsl:call-template name="getFlowcellMismatchGraphsRawLocalPath"/>
</xsl:template>


<xsl:template name="getFlowcellMismatchCurvesPfLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="'all'"/>
    <xsl:param name="sampleId" select="'all'"/>
    <xsl:param name="barcodeId" select="'all'"/>
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/PfCycleMismatchFragments.html')"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchCurvesRawLocalPath">
    <xsl:param name="flowcellId" select="../../../@flowcell-id"/>
    <xsl:param name="projectId" select="'all'"/>
    <xsl:param name="sampleId" select="'all'"/>
    <xsl:param name="barcodeId" select="'all'"/>
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/RawCycleMismatchFragments.html')"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchCurvesPfGlobalPath">
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/')"/><xsl:call-template name="getFlowcellMismatchCurvesPfLocalPath"/>
</xsl:template>

<xsl:template name="getFlowcellMismatchCurvesRawGlobalPath">
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_HTML_PARAM, '/')"/><xsl:call-template name="getFlowcellMismatchCurvesRawLocalPath"/>
</xsl:template>

<!-- Context: /Stats/Flowcell/Lane/Tile -->
<xsl:template name="getFlowcellLaneTileImageLocalPath">
    <xsl:param name="flowcellId" select="../../@flowcell-id"/>
    <xsl:param name="projectId" select="'all'"/>
    <xsl:param name="sampleId" select="'all'"/>
    <xsl:param name="barcodeId" select="'all'"/>
    <xsl:param name="laneNumber" select="../@number"/>
    <xsl:param name="tileNumber" select="@number"/>
    <xsl:param name="imageFileSuffix"/>
            
    <xsl:variable name="paddedTile" select="format-number($tileNumber, '0000')"/>
    
    <xsl:value-of select="concat($flowcellId, '/', $projectId, '/', $sampleId, '/', $barcodeId, '/s_', $laneNumber, '_', $paddedTile, '_', $imageFileSuffix)"/>
</xsl:template>

<xsl:template name="getFlowcellLaneTileImageGlobalPath">
    <xsl:param name="imageFileSuffix"/>

    <xsl:variable name="localPath">
        <xsl:call-template name="getFlowcellLaneTileImageLocalPath">
            <xsl:with-param name="imageFileSuffix" select="$imageFileSuffix"/>
        </xsl:call-template>
    </xsl:variable>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_IMAGES_PARAM, '/', $localPath)"/>
</xsl:template>


<xsl:template name="getFlowcellLaneTileMismatchesImageUniquePath">
    <xsl:param name="pf"/>
    <xsl:call-template name="getFlowcellLaneTileImageLocalPath">
        <xsl:with-param name="imageFileSuffix" select="concat($pf, $MISMATCHES_FILE_SUFFIX)"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchesImageFilePath">
    <xsl:param name="pf"/>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_IMAGES_PARAM, '/')"/><xsl:call-template name="getFlowcellLaneTileMismatchesImageUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchesThumbnailUniquePath">
    <xsl:param name="pf"/>
    <xsl:call-template name="getFlowcellLaneTileImageLocalPath">
        <xsl:with-param name="imageFileSuffix" select="concat($pf, $MISMATCHES_THUMBNAIL_FILE_SUFFIX)"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchesThumbnailFilePath">
    <xsl:param name="pf"/>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_IMAGES_PARAM, '/')"/><xsl:call-template name="getFlowcellLaneTileMismatchesThumbnailUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template>
</xsl:template>


<xsl:template name="getFlowcellLaneTileMismatchCurvesImageUniquePath">
    <xsl:param name="pf"/>
    <xsl:call-template name="getFlowcellLaneTileImageLocalPath">
        <xsl:with-param name="imageFileSuffix" select="concat($pf, $MISMATCH_CURVES_FILE_SUFFIX)"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchCurvesImageFilePath">
    <xsl:param name="pf"/>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_IMAGES_PARAM, '/')"/><xsl:call-template name="getFlowcellLaneTileMismatchCurvesImageUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchCurvesThumbnailUniquePath">
    <xsl:param name="pf"/>
    <xsl:call-template name="getFlowcellLaneTileImageLocalPath">
        <xsl:with-param name="imageFileSuffix" select="concat($pf, $MISMATCH_CURVES_THUMBNAIL_FILE_SUFFIX)"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="getFlowcellLaneTileMismatchCurvesThumbnailFilePath">
    <xsl:param name="pf"/>
    <xsl:value-of select="concat($OUTPUT_DIRECTORY_IMAGES_PARAM, '/')"/><xsl:call-template name="getFlowcellLaneTileMismatchCurvesThumbnailUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template>
</xsl:template>

</xsl:stylesheet>
