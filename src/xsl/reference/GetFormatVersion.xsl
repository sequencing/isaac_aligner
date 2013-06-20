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
 ** <https://github.com/sequencing/licenses/>.
 **
 ** The distribution includes the code libraries listed below in the
 ** 'redist' sub-directory. These are distributed according to the
 ** licensing terms governing each library.
 **
 ** \file GetFormatVersion.xsl
 **
 ** \author Roman Petrovski
 **/
-->
<xsl:stylesheet version="1.0" 
xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
>

<xsl:output method="text"/>

<xsl:template match="/"> 
    <xsl:value-of select="SortedReference/FormatVersion"/>
</xsl:template>

</xsl:stylesheet>
