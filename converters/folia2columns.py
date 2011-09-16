#!/usr/bin/env python
#-*- coding:utf-8 -*-

import getopt
import codecs
import sys
try:
    from pynlpl.formats import folia
except:
    print >>sys.stderr,"ERROR: pynlpl not found, please obtain PyNLPL from github: $ git clone git://github.com/proycon/pynlpl.git"

def usage():
    print >>sys.stderr, "folia2columns.py"
    print >>sys.stderr, "  by Maarten van Gompel (proycon)"
    print >>sys.stderr, "  Tilburg University / Radboud University Nijmegen"
    print >>sys.stderr, "  2011 - Licensed under GPLv3"
    print >>sys.stderr, ""
    print >>sys.stderr, "This conversion script reads a FoLiA XML document and produces a"
    print >>sys.stderr, "simple columned output format in which each token appears on one"
    print >>sys.stderr, "line. Note that only simple token annotations are supported and a lot"
    print >>sys.stderr, "of FoLiA data can not be intuitively expressed in a simple columned format!"
    print >>sys.stderr, ""
    print >>sys.stderr, "Parameters:"
    print >>sys.stderr, "  -f [FoLiA XML file]          Specify a FoLiA document to process (mandatory)"
    
    print >>sys.stderr, "  -c [columns]                 Comma separated list of desired column layout (mandatory), choose from:"
    print >>sys.stderr, "                               id      - output word ID"
    print >>sys.stderr, "                               text    - output the text of the word (the word itself)"    
    print >>sys.stderr, "                               pos     - output PoS annotation class"
    print >>sys.stderr, "                               poshead - output PoS annotation head feature"
    print >>sys.stderr, "                               lemma   - output lemma annotation class"
    print >>sys.stderr, "                               sense   - output sense annotation class"
    print >>sys.stderr, "                               phon    - output phonetic annotation class"
    print >>sys.stderr, "                               senid   - output sentence ID"
    print >>sys.stderr, "                               parid   - output paragraph ID"
    print >>sys.stderr, "                               N     - word/token number (absolute)"
    print >>sys.stderr, "                               n     - word/token number (relative to sentence)"
    print >>sys.stderr, "Options:"
    print >>sys.stderr, "  --csv                        Output in CSV format"            
    print >>sys.stderr, "  -o [filename]                Output to file instead of stdout"
    print >>sys.stderr, "  -e [encoding]                Output encoding (default: utf-8)"
    print >>sys.stderr, "  -H                           Suppress header output"    
    print >>sys.stderr, "  -S                           Suppress sentence spacing  (no whitespace between sentences)"    
    print >>sys.stderr, "  --csv                        Output in CSV format"        
    print >>sys.stderr, "  -x [sizeinchars]             Space columns for human readability (instead of plain tab-separated columns)"
    
try:
    opts, args = getopt.getopt(sys.argv[1:], "f:o:hHSc:x:", ["help", "csv"])
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(2)

filename = None

output_header = True
csv = False
outputfile = None
sentencespacing = True
nicespacing = 0

encoding = 'utf-8'
columnconf = []

for o, a in opts:
    if o == '-f':
        filename = a
    elif o == '-c':
        for a in a.split(','):
            columnconf.append(a)
    elif o == '-h':
        usage()
        sys.exit(0)
    elif o == '-H':
        output_header = False        
    elif o == '-S':
        sentencespacing = False
    elif o == '-e':
        encoding = a
    elif o == '-o':
        outputfile = a
    elif o == '-x':        
        nicespacing = int(a)
    elif o == '--csv':
        csv = True
    else:
        raise Exception("No such option: " + o)

if not filename:
    print >>sys.stderr,"ERROR: No FoLiA document specified (use -f)"
    usage()
    sys.exit(2)
elif not columnconf:
    print >>sys.stderr,"ERROR: No column configuration specified (use -c)"
    usage()
    sys.exit(2)
    
doc = folia.Document(file=filename)

prevsen = None


if outputfile: outputfile = codecs.open(outputfile,'w',encoding)
    
    
    

    
if nicespacing:

    def resize(s, i, spacing):
        if len(s) >= spacing[i]:
            s = s[0:spacing[i] - 1] + ' '
        elif len(s) < spacing[i]:
            s = s + (' ' * (spacing[i] - len(s)))
        #print '[' + s + ']', len(s), spacing[i]
        return s
    
    spacing = []
    for c in columnconf:
        if c == 'n':
            spacing.append(3)
        elif c == 'N':
            spacing.append(7)
        elif c == 'poshead':
            spacing.append(5)
        else:
            spacing.append(nicespacing)
    
if output_header:        
    
    if csv:
        columns = [ '"' + x.upper()  + '"' for x in columnconf ]
    else:
        columns = [ x.upper()  for x in columnconf ]

    if nicespacing and not csv:
        columns = [ resize(x, i, spacing) for i, x in enumerate(columnconf) ]
    
    if csv:
        line = ','.join(columns)
    else:
        line = '\t'.join(columns)

    if outputfile:
        outputfile.write(line + '\n')
    else:    
        print line.encode('utf-8')    

wordnum = 0



for i, w in enumerate(doc.words()):
    if w.sentence() != prevsen and i > 0:
        if sentencespacing:
            if outputfile:
                outputfile.write('\n')
            else:
                print         
        wordnum = 0
    prevsen = w.sentence()
    wordnum += 1
    columns = []
    for c in columnconf:
        if c == 'id':
            columns.append(w.id)
        elif c == 'text':
            columns.append(w.text())
        elif c == 'n':
            columns.append(str(wordnum))
        elif c == 'N':
            columns.append(str(i+1))
        elif c == 'pos':
            try:
                columns.append(w.annotation(folia.PosAnnotation).cls)    
            except:
                columns.append('-')
        elif c == 'poshead':
            try:
                columns.append(w.annotation(folia.PosAnnotation).feat('head'))    
            except:
                columns.append('-')     
        elif c == 'lemma':
            try:
                columns.append(w.annotation(folia.LemmaAnnotation).cls)    
            except:
                columns.append('-')
        elif c == 'sense':
            try:
                columns.append(w.annotation(folia.SenseAnnotation).cls)    
            except:
                columns.append('-')
        elif c == 'phon':
            try:
                columns.append(w.annotation(folia.PhonAnnotation).cls)    
            except:
                columns.append('-')                
        elif c == 'senid':
            columns.append(w.sentence().id)
        elif c == 'parid':            
            try:
                columns.append(w.paragraph().id)
            except:
                columns.append('-')
        elif c:
            print >>sys.stderr,"ERROR: Unsupported configuration: " + c
            sys.exit(1)
        
    if nicespacing and not csv:
        columns = [ resize(x,j, spacing) for j,x  in enumerate(columns) ]
        
    if csv:
        line = ",".join([ '"' + x  + '"' for x in columns ]) 
    else:
        line = "\t".join(columns)
        
    if outputfile:
        outputfile.write(line+'\n')
    else:    
        print line.encode('utf-8')

if outputfile:
    outputfile.close()
