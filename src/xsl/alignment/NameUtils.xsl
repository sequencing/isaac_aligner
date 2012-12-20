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
 ** \file NameUtils.xsl
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
<xsl:template name="getNodeDisplayName">
    <xsl:param name="node" select="."/>
    
    <xsl:choose>
        <xsl:when test="name()='Flowcell'">
            <xsl:choose>
                <xsl:when test="@flowcell-id='all'">[all flowcells]</xsl:when>
                <xsl:otherwise><xsl:value-of select="@flowcell-id"/></xsl:otherwise>
            </xsl:choose>
        </xsl:when>        <xsl:when test="name()='Project'">
            <xsl:choose>
                <xsl:when test="@name='all'">[all projects]</xsl:when>
                <xsl:when test="@name='default'">[default project]</xsl:when>
                <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
            </xsl:choose>
        </xsl:when>
        <xsl:when test="name()='Sample'">
            <xsl:choose>
                <xsl:when test="@name='all'">[all samples]</xsl:when>
                <xsl:when test="@name='unknown'">[unknown sample]</xsl:when>
                <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
            </xsl:choose>
        </xsl:when>
        <xsl:when test="name()='Barcode'">
            <xsl:choose>
                <xsl:when test="@name='all'">[all barcodes]</xsl:when>
                <xsl:when test="@name='unknown'">[unknown barcode]</xsl:when>
                <xsl:otherwise><xsl:value-of select="@name"/></xsl:otherwise>
            </xsl:choose>
        </xsl:when>

        <xsl:otherwise>
                <xsl:value-of select="@name"/>
        </xsl:otherwise>
    </xsl:choose>
</xsl:template>

<xsl:template name="getBarcodeDisplayPath">
    <p>
        <xsl:for-each select="../../.."><xsl:call-template name="getNodeDisplayName"/></xsl:for-each> /
        <xsl:for-each select="../.."><xsl:call-template name="getNodeDisplayName"/></xsl:for-each> /
        <xsl:for-each select=".."><xsl:call-template name="getNodeDisplayName"/></xsl:for-each> /
        <xsl:for-each select="."><xsl:call-template name="getNodeDisplayName"/></xsl:for-each>
    </p>
</xsl:template>

</xsl:stylesheet>
