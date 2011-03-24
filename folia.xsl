<?xml version="1.0" encoding="utf-8" ?>
<xsl:stylesheet version="1.0" xmlns="http://www.w3.org/1999/xhtml" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:imdi="http://www.mpi.nl/IMDI/Schema/IMDI" xmlns:folia="http://ilk.uvt.nl/folia">

<xsl:output method="html" encoding="UTF-8" omit-xml-declaration="yes" doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN" doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd" indent="yes" cdata-section-elements="script"/>


<xsl:template match="/folia:FoLiA">
    <html>
        <head>
            <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
            <link rel="stylesheet" href="http://ilk.uvt.nl/folia/folia.css" type="text/css"></link>
            <script type="text/javascript" src="http://code.jquery.com/jquery-1.5.1.min.js"></script>
            <xsl:choose>
             <xsl:when test="/FoLiA/metadata/imdi:METATRANSCRIPT/imdi:Session/imdi:Title">
                <title><xsl:value-of select="/FoLiA/metadata/imdi:METATRANSCRIPT/imdi:Session/imdi:Title"/></title>
             </xsl:when>
             <xsl:otherwise>
                <title><xsl:value-of select="@xml:id"/></title>
             </xsl:otherwise>
            </xsl:choose>            
        </head>
        <body>
            <xsl:apply-templates select="//imdi:METATRANSCRIPT" />
            <xsl:apply-templates select="folia:text" />
        </body>
    </html>
</xsl:template>

<xsl:template match="//imdi:METATRANSCRIPT">
    <div id="metadata">
        <table>
        <tr><th>Name:</th><td><xsl:value-of select="imdi:Session/imdi:Name"/></td></tr>
        <tr><th>Title:</th><td><strong><xsl:value-of select="imdi:Session/imdi:Title"/></strong></td></tr>
        <tr><th>Date:</th><td><xsl:value-of select="imdi:Session/imdi:Date"/></td></tr>
        <xsl:if test="//imdi:Source/imdi:Access/imdi:Availability">
                <tr><th>Availability:</th><td><xsl:value-of select="//imdi:Source/imdi:Access/imdi:Availability"/></td></tr>
        </xsl:if>
        <xsl:if test="//imdi:Source/imdi:Access/imdi:Publisher">
                <tr><th>Publisher:</th><td><xsl:value-of select="//imdi:Source/imdi:Access/imdi:Publisher"/></td></tr>
        </xsl:if>
        </table>
    </div>
</xsl:template>


<xsl:template match="folia:text">
 <div class="text">
   <xsl:choose>
   <xsl:when test="/folia:div">
    <xsl:apply-templates select="/folia:div" />
   </xsl:when>
   <xsl:when test="//folia:p">
    <xsl:apply-templates select="//folia:p|//folia:head" />
   </xsl:when>
   <xsl:when test="//folia:s">
    <xsl:apply-templates select="//folia:s|//folia:head" />
   </xsl:when> 
   <xsl:otherwise>
    <span class="error">No content found in this text!</span>
   </xsl:otherwise>
  </xsl:choose>
 </div>
</xsl:template>

<xsl:template match="folia:div">
 <div class="div">
  <xsl:apply-templates />
 </div>
</xsl:template>

<xsl:template match="folia:p">
 <p>
  <xsl:apply-templates />
 </p>
</xsl:template>


<xsl:template match="folia:head">
 <h1>
  <xsl:apply-templates />
 </h1>
</xsl:template>

<xsl:template match="folia:s">
 <span class="s">
  <xsl:apply-templates />
 </span>
</xsl:template>

<xsl:template match="folia:w">
 <span id="{@xml:id}" class="word">
        <span class="attributes">
                <span class="wordid"><xsl:value-of select="@xml:id" /></span>
                <dl>
                        <xsl:apply-templates />
                </dl>
        </span>
        <xsl:value-of select="folia:t"/>
 </span>
 <xsl:text> </xsl:text> <!-- TODO: implement @nospace check -->
</xsl:template>

<xsl:template match="folia:pos">
 <dt>PoS</dt><dd><xsl:value-of select="@folia:class"/></dd>
</xsl:template>

<xsl:template match="folia:lemma">
 <dt>Lemma</dt><dd><xsl:value-of select="@folia:class"/></dd>
</xsl:template>

</xsl:stylesheet>
