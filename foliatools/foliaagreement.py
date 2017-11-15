#!/usr/bin/env python
#-*- coding:utf-8 -*-

from __future__ import print_function, unicode_literals, division, absolute_import

import argparse
import sys
from collections import Counter, defaultdict
try:
    from pynlpl.formats import folia
except:
    print("ERROR: pynlpl not found, please obtain PyNLPL from the Python Package Manager ($ sudo pip install pynlpl) or directly from github: $ git clone git://github.com/proycon/pynlpl.git",file=sys.stderr)
    sys.exit(2)

def get_corrections(doc, Class, foliaset):
    """Get relevant corrections from document"""
    for correction in doc.select(folia.Correction):
        if correction.hasnew():
            for annotation in correction.new():
                structural = isinstance(annotation, folia.AbstractStructure) #structural correction
                target = None
                if structural:
                    if annotation.hasannotation(Class, foliaset):
                        target = annotation
                    #TODO: deal with deletions
                elif isinstance(annotation, Class) and annotation.set == foliaset: #normal correction
                    if isinstance(annotation, folia.AbstractSpanAnnotation):
                        target = tuple(annotation.wrefs()) #TODO: implement handling
                    else:
                        target = annotation.ancestor(folia.AbstractStructure)
                if target is not None:
                    yield correction, structural, annotation.annotation(Class, foliaset), target

def inter_annotator_agreement(docs, Class, foliaset, corrections=False, verbose=False):
    nr = len(docs)
    index = []
    if corrections:
        for i, doc in enumerate(docs):
            index.append(defaultdict(dict))
            for correction, structural, annotation, target in get_corrections(doc, Class, foliaset):
                if target in index[i]:
                    print("WARNING: Overlapping annotation for " + repr(target) + ", overwriting!",file=sys.stderr)
                index[i][target] = annotation
    else:
        raise NotImplementedError #TODO

    #linking step: links annotations on the same things
    links = []
    targets = []
    for target in index[0]:
        linkchain = []
        for i in range(0,len(index)):
            if i == 0: targets.append(target)
            if target not in index[i]:
                break
            else:
                linkchain.append(index[i][target])
        if len(linkchain) == nr:
            links.append(linkchain)


    #evaluation step
    weakmatches = len(links) #match regardless of class

    strongmatches = 0 #match including class
    #compute strong matches
    for target, linkchain in zip(targets, links):
        match = all_equal([ get_value(annotation, Class) for annotation in linkchain ])
        strongmatches += int(match)
        if verbose:
            if match:
                print("YES\t" + target.id + "\t" + get_value(annotation, Class))
            else:
                print("NO\t" + target.id + "\t" + get_value(annotation, Class))

    #collect all possible targets (for normalisation)
    alltargets = set()
    for i in index:
        for target in index[i]:
            alltargets.add(hash(target))

    return strongmatches, weakmatches, len(alltargets)

def all_equal(collection):
    iterator = iter(collection)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)

def get_value(annotation, Class):
    if Class is folia.TextContent or Class is folia.PhonContent:
        return str(annotation)
    return annotation.cls

def main():
    parser = argparse.ArgumentParser(description="FoLiA Inter-Annotator Agreement: This tool computes inter-annotator agreement on two or more structurally equivalent FoLiA documents. Evaluation is expressed as accuracy on the total number of annotation targets (often words) and comes in two flavours: weak and strong. Weak checks only if the same items were marked and can be used as a measure of detection; strong checks if the assigned classes are equal amongst annotators.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.add_argument('--storeconst',dest='settype',help="", action='store_const',const='somevalue')
    parser.add_argument('-t','--type', type=str,help="Annotation type to consider", action='store',default="",required=True)
    parser.add_argument('-s','--set', type=str,help="Set definition (required if there is ambiguity in the document)", action='store',default="",required=False)
    parser.add_argument('-c','--corrections', help="Use corrections", action='store_true',default="",required=False)
    parser.add_argument('-v','--verbose', help="Verbose, list all matches/mismatches", action='store_true',required=False)
    #parser.add_argument('-i','--number',dest="num", type=int,help="", action='store',default="",required=False)
    parser.add_argument('documents', nargs='+', help='FoLiA Documents')
    args = parser.parse_args()

    for docfile in args.documents:
        docs = folia.Document(file=docfile)

    try:
        Type = folia.XML2CLASS[args.type]
    except KeyError:
        print("No such type: ", args.type,file=sys.stderr)

    foliaset = args.set

    strongmatches, weakmatches, total = inter_annotator_agreement(docs, Type, foliaset, args.corrections, args.verbose)
    if not total:
        print("strong\t0\t0")
        print("weak\t0\t0")
        print("total\t0")
    else:
        print("strong\t" + str(strongmatches) + "\t" + str(round(strongmatches/total,3)))
        print("weak\t" + str(weakmatches) + "\t" + str(round(weakmatches/total,3)))
        print("total\t" + str(total))

if __name__ == "__main__":
    main()
