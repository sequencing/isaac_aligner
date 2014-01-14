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
 ** \file TileMismatchPages.xsl
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

<xsl:key name="tilesById" match="/Stats/Flowcell/Lane/Tile" use="@number"/>

<xsl:template name="generateTileTable">
    <xsl:param name="flowcellNode"/>
    <xsl:param name="pf"/>
    <xsl:param name="imageFileSuffix"/>
    <xsl:param name="thumbnailFileSuffix"/>
    
<table border="1" cellpadding="5">
        <th>Tile:</th>
        <xsl:for-each select="$flowcellNode/Lane">
        <!-- xsl:for-each select="$flowcellNode/Lane/Tile"-->
        <!-- <xsl:for-each select="$flowcellNode/Lane/Tile[(key('tilesById', @number)[../../@flowcell-id=$flowcellNode/@flowcell-id])[1]/../@number=../@number]"> -->
        <xsl:sort select="@number" data-type="number"/>
            <th>
                <xsl:value-of select="concat('lane', @number)"/>
            </th>
        </xsl:for-each>
    
    <!-- <xsl:for-each select="Tile"> -->
    <xsl:for-each select="$flowcellNode/Lane/Tile[(key('tilesById', @number)[../../@flowcell-id=$flowcellNode/@flowcell-id])[1]/../@number=../@number]">
        <xsl:sort select="@number" data-type="number"/>
        <xsl:variable name="tileNumber" select="@number"/>

        <tr>
        <td><xsl:value-of select="$tileNumber"/></td>
        
        <xsl:for-each select="$flowcellNode/Lane">
            <xsl:sort select="@number" data-type="number"/>
            <xsl:variable name="laneNumber" select="@number"/>
            
            <xsl:variable name="thumbnailFilePath">
                <xsl:call-template name="getFlowcellLaneTileImageGlobalPath">
                    <xsl:with-param name="imageFileSuffix" select="concat($pf, $thumbnailFileSuffix)"/>
                    <xsl:with-param name="flowcellId" select="$flowcellNode/@flowcell-id"/>
                    <xsl:with-param name="laneNumber" select="$laneNumber"/>
                    <xsl:with-param name="tileNumber" select="$tileNumber"/>
                </xsl:call-template>
            </xsl:variable>
            <xsl:variable name="imageFilePath">
                <xsl:call-template name="getFlowcellLaneTileImageGlobalPath">
                    <xsl:with-param name="flowcellId" select="$flowcellNode/@flowcell-id"/>
                    <xsl:with-param name="laneNumber" select="$laneNumber"/>
                    <xsl:with-param name="tileNumber" select="$tileNumber"/>
                    <xsl:with-param name="imageFileSuffix" select="concat($pf, $imageFileSuffix)"/>
                </xsl:call-template>
            </xsl:variable>
            
            <td>
                <xsl:element name="a">
                    <xsl:attribute name="href"><xsl:value-of select="$imageFilePath"/></xsl:attribute>
                    <xsl:element name="img">
                        <xsl:attribute name="height">84</xsl:attribute>
                        <xsl:attribute name="width">84</xsl:attribute>
                        <xsl:attribute name="src"><xsl:value-of select="$thumbnailFilePath"/></xsl:attribute>
                    </xsl:element>
                </xsl:element>
            </td>
        </xsl:for-each>
        </tr>
    </xsl:for-each>
</table>
    
</xsl:template>

<xsl:template name="generateTileMismatchGraphsGnuplotScript">
    <xsl:param name="flowcellNode"/>
    <xsl:param name="pf"/>

    
    <xsl:for-each select="$flowcellNode/Lane">
    <xsl:sort select="@number"/>
        <xsl:variable name="lane" select="@number"/>
        <xsl:for-each select="Tile">
        <xsl:sort select="@number"/>
            <xsl:variable name="thumbnailFileName"><xsl:call-template name="getFlowcellLaneTileMismatchesThumbnailUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="thumbnailFilePath"><xsl:call-template name="getFlowcellLaneTileMismatchesThumbnailFilePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="imageFileName"><xsl:call-template name="getFlowcellLaneTileMismatchesImageUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="imageFilePath"><xsl:call-template name="getFlowcellLaneTileMismatchesImageFilePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            
<!-- image generation -->
            <xsl:variable name="maxX" select="math:max(*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle/@number)"/>
            <xsl:variable name="minX" select="math:min(*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle/@number)"/>
            <xsl:variable name="imageWidth" select="($maxX - $minX + 1) * 3"/>
            <xsl:variable name="lineWidth" select="0.5"/><!--fixed($maxX div $imageWidth)"/-->
reset
set terminal gif medium size <xsl:value-of select="$imageWidth"/>,600 \
    &#35;ffffff &#35;000000 &#35;404040 &#35;ff0000 &#35;00ff00 &#35;0000ff &#35;000000 &#35;000000 &#35;0000ff &#35;000000 &#35;000000 &#35;777777;  
set output '<xsl:value-of select="$imageFilePath"/>';
set x2range [<xsl:value-of select="$minX"/>:<xsl:value-of select="$maxX"/>];
set xrange [<xsl:value-of select="$minX"/>:<xsl:value-of select="$maxX"/>];
set yrange [0:20];
set ytics nomirror 5;
unset xlabel;
set x2label "Cycle Number";
unset xtics;
set x2tics;
set grid ytics x2tics
set boxwidth <xsl:value-of select="$lineWidth"/>;

set title "<xsl:value-of select="$imageFileName"/> \n Uniquely aligned fragments: \
<xsl:for-each select="*[name()=$pf]/Read">
<xsl:value-of select="concat('R', @number, ':', format-number(UniquelyAlignedFragments/Count, '###,###,###,###,###'), ' ')"/>\
</xsl:for-each>";
set nokey;
set style fill solid 1.0 noborder;
set size 1.0, 1.0;

set multiplot
set size 1.0, 0.5;
unset xlabel;
set origin 0.0, 0.5;
set bmargin 0;
unset tmargin;


set ylabel "% mismatches";
plot '-' using 1:2 with boxes lt 1;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', format-number(100 * (@mismatches) div $uniqAlignFragments, '#.#'))"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>

set yrange [20:0];
unset x2label;
set xlabel "Cycle Number";
set xtics;
unset x2tics;
set grid ytics xtics

set origin 0.0, 0.0;
set tmargin 0;
unset bmargin;

set ylabel "% blanks";
unset title
plot '-' using 1:2 with boxes lt 3;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', format-number(100 * (@blanks) div $uniqAlignFragments, '#.#'))"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
unset multiplot

<!-- thumbnail generation -->
reset
set terminal gif medium size 84,84 \
    &#35;ffffff &#35;000000 &#35;404040 &#35;ff0000 &#35;00ff00 &#35;0000ff &#35;000000 &#35;000000 &#35;0000ff &#35;000000 &#35;000000 &#35;777777;  
set output '<xsl:value-of select="$thumbnailFilePath"/>';
set xrange [<xsl:value-of select="$minX"/>:<xsl:value-of select="$maxX"/>];
set yrange [0:20];
unset ytics;
unset xlabel;
unset xtics;
set nokey;
set size 1.0, 1.0;
set lmargin 0;
set rmargin 0;
set tmargin 0;
set bmargin 0;

set multiplot
set size 1.0, 0.5;
set origin 0.0, 0.5;

plot '-' using 1:2 with impulses lt 1;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', format-number(100 * (@mismatches) div $uniqAlignFragments, '#.#'))"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>

set yrange [20:0];
set origin 0.0, 0.0;
plot '-' using 1:2 with impulses lt 3;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/MismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', format-number(100 * (@blanks) div $uniqAlignFragments, '#.#'))"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
unset multiplot

        </xsl:for-each>
    </xsl:for-each>

</xsl:template>


<xsl:template name="generateTileMismatchGraphsPage">
    <xsl:param name="flowcellNode"/>
    <xsl:param name="pf"/>
    
<html>
<link rel="stylesheet" href="../../../../{$CSS_FILE_NAME}" type="text/css"/>
<body>

    <xsl:call-template name="generateTileTable">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="$pf"/>
        <xsl:with-param name="imageFileSuffix" select="$MISMATCHES_FILE_SUFFIX"/>
        <xsl:with-param name="thumbnailFileSuffix" select="$MISMATCHES_THUMBNAIL_FILE_SUFFIX"/>
    </xsl:call-template>

<p>@iSAAC_VERSION_FULL@</p>
</body>
</html>

</xsl:template>


<xsl:template name="generateTileMismatchGraphsPfPage">
    <xsl:param name="flowcellNode"/>

    <xsl:call-template name="generateTileMismatchGraphsPage">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="'Pf'"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="generateTileMismatchGraphsRawPage">
    <xsl:param name="flowcellNode"/>

    <xsl:call-template name="generateTileMismatchGraphsPage">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="'Raw'"/>
    </xsl:call-template>
</xsl:template>


<xsl:template name="generateTileMismatchCurvesGnuplotScript">
    <xsl:param name="flowcellNode"/>
    <xsl:param name="pf"/>

    
    <xsl:for-each select="$flowcellNode/Lane">
    <xsl:sort select="@number"/>
        <xsl:variable name="lane" select="@number"/>
        <xsl:for-each select="Tile">
        <xsl:sort select="@number"/>
        
            <xsl:variable name="thumbnailFileName"><xsl:call-template name="getFlowcellLaneTileMismatchCurvesThumbnailUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="thumbnailFilePath"><xsl:call-template name="getFlowcellLaneTileMismatchCurvesThumbnailFilePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="imageFileName"><xsl:call-template name="getFlowcellLaneTileMismatchCurvesImageUniquePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
            <xsl:variable name="imageFilePath"><xsl:call-template name="getFlowcellLaneTileMismatchCurvesImageFilePath"><xsl:with-param name="pf" select="$pf"/></xsl:call-template></xsl:variable>
        
<!-- image generation -->
            <xsl:variable name="maxX" select="math:max(*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle/@number)"/>
            <xsl:variable name="minX" select="math:min(*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle/@number)"/>
            <xsl:variable name="imageWidth" select="($maxX - $minX + 1) * 3"/>
            <xsl:variable name="lineWidth" select="0.5"/><!--fixed($maxX div $imageWidth)"/-->
reset
set terminal gif medium size <xsl:value-of select="$imageWidth"/>,600 \
    &#35;ffffff &#35;000000 &#35;404040 &#35;ff0000 &#35;00ff00 &#35;0000ff &#35;000000 &#35;000000 &#35;0000ff &#35;000000 &#35;000000 &#35;777777;  
set output '<xsl:value-of select="$imageFilePath"/>';
set xrange [<xsl:value-of select="$minX"/>:<xsl:value-of select="$maxX"/>];
set yrange [0:100];
set ytics nomirror 10;
set xtics nomirror;
set boxwidth <xsl:value-of select="$lineWidth"/>;

set title "<xsl:value-of select="$imageFileName"/> \n Uniquely aligned fragments: \
<xsl:for-each select="*[name()=$pf]/Read">
<xsl:value-of select="concat('R', @number, ':', format-number(UniquelyAlignedFragments/Count, '###,###,###,###,###'), ' ')"/>\
</xsl:for-each>";
set xlabel "Cycle Number";
set ylabel "% uniquely aligned fragments with 'x' mismatches or less";
set key left bottom below
set style fill solid 1.0 noborder

plot '-' using 1:2 title '4 or less' with boxes lt 7, \
     '-' using 1:2 title '3 or less' with boxes lt 2, \
     '-' using 1:2 title '2 or less' with boxes lt 3, \
     '-' using 1:2 title '1 or less' with boxes lt 1, \
     '-' using 1:2 title '0 mismatches' with boxes lt 9;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @four) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @three) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @two) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @one) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>

<!-- thumbnail generation -->
reset
set terminal gif medium size 84,84 \
    &#35;ffffff &#35;000000 &#35;404040 &#35;ff0000 &#35;00ff00 &#35;0000ff &#35;000000 &#35;000000 &#35;0000ff &#35;000000 &#35;000000 &#35;777777;  
set output '<xsl:value-of select="$thumbnailFilePath"/>';
set xrange [<xsl:value-of select="$minX"/>:<xsl:value-of select="$maxX"/>];
set yrange [0:100];
set noytics;
set noxtics;
set nokey;
set notitle;
set noxlabel;
set noylabel;
set lmargin 0;
set rmargin 0;
set tmargin 0;
set bmargin 0;

plot '-' using 1:2 with impulses lt 7, \
     '-' using 1:2 with impulses lt 2, \
     '-' using 1:2 with impulses lt 3, \
     '-' using 1:2 with impulses lt 1, \
     '-' using 1:2 with impulses lt 9;
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @four) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @three) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @two) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more + @one) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>
<xsl:for-each select="*[name()=$pf]/Read/UniquelyAlignedFragments/FragmentMismatchesByCycle/Cycle">
<xsl:sort select="@number" type="number"/>
<xsl:variable name="uniqAlignFragments" select="../../Count"/>
<xsl:value-of select="concat(@number, ' ', 100 * ($uniqAlignFragments - @more) div $uniqAlignFragments)"/>
<xsl:text>
</xsl:text>
</xsl:for-each><xsl:text>e
</xsl:text>

        </xsl:for-each>
    </xsl:for-each>

</xsl:template>


<xsl:template name="generateTileMismatchCurvesPage">
    <xsl:param name="flowcellNode"/>
    <xsl:param name="pf"/>

<html>
<link rel="stylesheet" href="../{$CSS_FILE_NAME}" type="text/css"/>
<body>

    <xsl:call-template name="generateTileTable">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="$pf"/>
        <xsl:with-param name="imageFileSuffix" select="$MISMATCH_CURVES_FILE_SUFFIX"/>
        <xsl:with-param name="thumbnailFileSuffix" select="$MISMATCH_CURVES_THUMBNAIL_FILE_SUFFIX"/>
    </xsl:call-template>

<p>@iSAAC_VERSION_FULL@</p>
</body>
</html>

</xsl:template>

<xsl:template name="generateTileMismatchCurvesPfPage">
    <xsl:param name="flowcellNode"/>

    <xsl:call-template name="generateTileMismatchCurvesPage">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="'Pf'"/>
    </xsl:call-template>
</xsl:template>

<xsl:template name="generateTileMismatchCurvesRawPage">
    <xsl:param name="flowcellNode"/>

    <xsl:call-template name="generateTileMismatchCurvesPage">
        <xsl:with-param name="flowcellNode" select="$flowcellNode"/>
        <xsl:with-param name="pf" select="'Raw'"/>
    </xsl:call-template>
</xsl:template>


</xsl:stylesheet>
