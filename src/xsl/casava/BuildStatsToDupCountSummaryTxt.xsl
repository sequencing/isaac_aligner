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
 ** \file BuildStatsToDupCountSummaryTxt.xsl
 **
 ** \brief Translation from BuildStats.xml to CASAVA dupCount.summary.txt
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
<xsl:text># ** CASAVA reference sequence duplicate summary **
#$ COLUMNS seq_name nonduplicate_pairs duplicate_pairs
</xsl:text>
</xsl:variable>


<xsl:template match="/">
    <xsl:value-of select="$readsIdxHeader"/>
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="/Stats">

    <xsl:for-each select="Project[@name=$ISAAC_PROJECT_NAME_PARAM]/Sample[@name=$ISAAC_SAMPLE_NAME_PARAM]/Contig">
    <xsl:sort select='ReferenceTotalBases' data-type='number' order='descending'/>
    
        <xsl:variable name="seq_name" select="@name"/>
        
        <xsl:variable name="nonduplicate_pairs" select="sum(Bin/UniqueFragments) div 2"/>
        <xsl:variable name="duplicate_pairs" select="(sum(Bin/TotalFragments) - sum(Bin/UniqueFragments)) div 2"/>

        <xsl:value-of select="$seq_name"/>
        <xsl:value-of select="$tab"/>
        <xsl:value-of select="$nonduplicate_pairs"/>
        <xsl:value-of select="$tab"/>
        <xsl:value-of select="$duplicate_pairs"/>
        
        <xsl:value-of select="$newline"/>
    </xsl:for-each>

</xsl:template>

</xsl:stylesheet>
