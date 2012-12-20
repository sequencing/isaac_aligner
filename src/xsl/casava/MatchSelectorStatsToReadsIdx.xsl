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
 ** \file MatchSelectorStatsToReadsIdx.xsl
 **
 ** \brief Translation from MatchSelectorStats.xml to CASAVA Reads.idx 
 ** 
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:str="http://exslt.org/strings"
> 

<xsl:output method="text"/>
<xsl:strip-space elements="*"/>

<xsl:variable name="newline"><xsl:text>
</xsl:text></xsl:variable>

<xsl:variable name="tab"><xsl:text>&#x09;</xsl:text></xsl:variable>

<xsl:variable name="readsIdxHeader">
<xsl:text>#\$ COLUMNS totalReads&#x09;totalBases&#x09;usedReads&#x09;usedBases&#x09;goodReads&#x09;</xsl:text>
<xsl:text>goodReadBases&#x09;failedTileFilt&#x09;failedFilt&#x09;failedQC&#x09;singleExclude&#x09;nonUniqueAlign&#x09;</xsl:text>
<xsl:text>goodPairs&#x09;mixedPairs&#x09;riboCount&#x09;mitoCount&#x09;splicedReads&#x09;exportFileDir&#x09;</xsl:text>
<xsl:text>exportSet&#x09;date&#x09;instrumentName&#x09;runID&#x09;laneNumber&#x09;barcode&#x09;highSDFragmentLength&#x09;</xsl:text>
<xsl:text>lowSDFragmentLength&#x09;maxFragmentLength&#x09;medianFragmentLength&#x09;minFragmentLength&#x09;</xsl:text>
<xsl:text>nominalOrientation&#x09;read1Length&#x09;read2Length
</xsl:text>
</xsl:variable>


<xsl:template match="/">
    <xsl:value-of select="$readsIdxHeader"/>
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="/Stats/Flowcell">

    <xsl:for-each select="Project[@name=$ISAAC_PROJECT_NAME_PARAM]/Sample[@name=$ISAAC_SAMPLE_NAME_PARAM]/Barcode[@name!='all']/Lane/Tile">
    
        <xsl:if test="'true' = ./AssumedTemplateLength/Stable">
            <xsl:variable name="totalReads" select="sum(./*/Read/AllFragments/Count)"/>
            <xsl:variable name="totalBases" select="sum(./*/Read/AllFragments/Yield)"/>
            <xsl:variable name="usedReads" select="sum(./Pf/Read/AllFragments/Count)"/>
            <xsl:variable name="usedBases" select="sum(./Pf/Read/AllFragments/Yield)"/>
            <xsl:variable name="goodReads" select="''"/><!-- depends on casava QVCutoffSingle-->
            <xsl:variable name="goodReadBases" select="''"/><!-- depends on casava QVCutoffSingle-->
            <xsl:variable name="failedTileFilt" select="0"/><!-- CASAVA tile filtering is deprecated in CASAVA-->
            <xsl:variable name="failedFilt" select="sum(./Raw/Read/AllFragments/Count) - sum(./Pf/Read/AllFragments/Count)"/>
            <xsl:variable name="failedQC" select="sum(./*/AlignmentScoreDistribution/Score[@score&lt;4]/@templates)"/>
            <xsl:variable name="unanchored" select="sum(./*/UnanchoredClusterCount)"/>
            <xsl:variable name="nonUniquePairs" select="''"/><!-- Currently iSAAC does not keep neighbourhood record for fragments-->
        
            <!-- in CASAVA depends QVCutoff, however being used in PairStats.cpp as 
            if it was the clusters to which the template length applies.
            iSAAC computes template length on PF clusters untill it settles. Just
            using the total Pf cluster count-->
            <xsl:variable name="goodPairs" select="sum(./Pf/ClusterCount)"/>
            
            <xsl:variable name="mixedPairs" select="''"/><!-- TODO: see if it is worth supporting-->
            <xsl:variable name="repeatMatchCount" select="''"/><!-- TODO: see if it is worth supporting-->
            <xsl:variable name="mitoCount" select="'0'"/><!-- RNA feature-->
            <xsl:variable name="splicedReads" select="''"/><!-- RNA feature-->
            <xsl:variable name="exportFileDir" select="''"/><!-- non-essential path to the source data-->
            <xsl:variable name="exportSet" select="@number"/><!-- pretend there is one tile per set-->
            <xsl:variable name="date" select="''"/><!-- non-essential ambiguous value-->
        
            <!-- TODO: get run id out of iSAAC parameters-->
            <xsl:variable name="instrumentName" select="substring-before(../../../../../@flowcell-id, '_')"/>
        
            <xsl:variable name="runID" select="substring-after(../../../../../@flowcell-id, '_')"/>
            <xsl:variable name="laneNumber" select="../@number"/>
            <xsl:variable name="barcodeSequence" select="../../@name"/>
            <xsl:variable name="highSDFragmentLength" select="./AssumedTemplateLength/HighStdDev"/>
            <xsl:variable name="lowSDFragmentLength" select="./AssumedTemplateLength/LowStdDev"/>
            <xsl:variable name="maxFragmentLength" select="./AssumedTemplateLength/Max"/>
            <xsl:variable name="medianFragmentLength" select="./AssumedTemplateLength/Median"/>
            <xsl:variable name="minFragmentLength" select="./AssumedTemplateLength/Min"/>
            <xsl:variable name="nominalOrientation" select="./AssumedTemplateLength/Class1"/>
            <xsl:variable name="read1SeqLength" select="../../../../../Read[@number=1]/Length"/>
            <xsl:variable name="read2SeqLength" select="../../../../../Read[@number=2]/Length"/>
                
            <xsl:value-of select="$totalReads"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$totalBases"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$usedReads"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$usedBases"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$goodReads"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$goodReadBases"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$failedTileFilt"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$failedFilt"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$failedQC"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$unanchored"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$nonUniquePairs"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$goodPairs"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$mixedPairs"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$repeatMatchCount"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$mitoCount"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$splicedReads"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$exportFileDir"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$exportSet"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$date"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$instrumentName"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$runID"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$laneNumber"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$barcodeSequence"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$highSDFragmentLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$lowSDFragmentLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$maxFragmentLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$medianFragmentLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$minFragmentLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$nominalOrientation"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$read1SeqLength"/>
            <xsl:value-of select="$tab"/>
            <xsl:value-of select="$read2SeqLength"/>
            
            <xsl:value-of select="$newline"/>
        </xsl:if>
    </xsl:for-each>
</xsl:template>

</xsl:stylesheet>
