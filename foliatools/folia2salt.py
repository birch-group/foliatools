#!/usr/bin/env python3

"""FoLiA to Salt XML conversion. Can in turn be used by Pepper for further conversion to other formats."""
import sys
import os
import argparse
import glob
from collections import OrderedDict
import lxml.etree
from lxml.builder import ElementMaker
from foliatools import VERSION as TOOLVERSION
import folia.main as folia
E = ElementMaker(nsmap={"sDocumentStructure":"sDocumentStructure", "xmi":"http://www.omg.org/XMI", "xsi": "http://www.w3.org/2001/XMLSchema-instance", "saltCore":"saltCore","saltCommon":"saltCommon", "sCorpusStructure":"sCorpusStructure" })
DESCRIPTION_TAGS = []

# Birch specific imports
from feats import FEAT_DICT
from Blogger import Blogger
import xmltodict as xd
logger = Blogger()

# Defining which nodes to ignore when adding layer info
BASE_LAYERS =  ['sDocumentStructure:STextualDS',
                'sDocumentStructure:SAudioDS',
                ]
"""
TODO:
    - Add support for corresponding audio files
        - From https://raw.githubusercontent.com/korpling/salt/master/gh-site/doc/salt_modelGuide.pdf:
        - 'To  create  an  audio  annotation, theSAnnotation.sValue must be set to an URI pointing to the audio file'
        - SAudioDSRelation behaves similar to STextualRelation tags
"""

#################################################
                # Tools #
#################################################
def etree_to_dict(t):
    d = {t.tag : map(etree_to_dict, t.iterchildren())}
    d.update(('@' + k, v) for k, v in t.attrib.iteritems())
    d['text'] = t.text
    return d

def processdir(d, **kwargs):
    print("Searching in  " + d,file=sys.stderr)
    corpusdocs = []
    extension = kwargs.get('extension','xml').strip('.')
    for ix, f in enumerate(glob.glob(os.path.join(d ,'*'))):
        if f[-len(extension) - 1:] == '.' + extension:
            _, corpusdoc = convert(f, **kwargs)
            corpusdocs.append(corpusdoc)
        elif kwargs.get('recurse') and os.path.isdir(f):
            corpusdocs += processdir(f,**kwargs)
    return corpusdocs

def time_tuple_to_string(t):
    """
    Formats Folia Doc time tuple to string for Salt Doc time relation.
    e.g. (0, 0, 2, 100) --> 00:00:02.100
    """
    out = ""
    for ix, item in enumerate(t):
        if ix == 3:
            out += f'.{str(item)}'
        elif len(str(item)) == 1:
            out += f':0{str(item)}'
        else:
            out += f':{str(item)}'
        if ix == 1:
            out = out[1:]
    return out

#################################################
                # Conversion Scripts #
#################################################
def convert(f, **kwargs):
    """Convert a FoLiA document to Salt XML"""
    success = False
    doc = folia.Document(file=f)
    if not doc.declared(folia.AnnotationType.TOKEN):
        raise Exception("Only tokenized documents can be handled at the moment")
    layers = OrderedDict()
    # return convert_tokens(doc, layers, **kwargs)
    nodes, textrelations, text, phonrelations, phon, map_id_to_nodenr, nodes_seqnr, node_ix = convert_tokens(doc, layers, **kwargs)

    edges = textrelations

    layers = list(build_layers(layers, nodes, edges, node_ix)) #this modifies the nodes and edges as well

    #Create the document (sDocumentGraph)
    saltdoc = getattr(E,"{sDocumentStructure}SDocumentGraph")(
        {"{http://www.omg.org/XMI}version":"2.0"},
        E.labels({ # document ID
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
            "namespace": "salt",
            "name": "id",
            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id
        }),
        *nodes,
        *edges,
        *layers)


    corpusdoc = E.nodes( #sDocument
                        {"{http://www.w3.org/2001/XMLSchema-instance}type": "sCorpusStructure:SDocument"},
                        E.labels({ # document ID
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                            "namespace": "salt",
                            "name": "id",
                            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SNAME",
                            "value": "T::" + doc.id
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SDOCUMENT_GRAPH",
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SDOCUMENT_GRAPH_LOCATION",
                            "value": "T::file:" + os.path.join(os.path.realpath(kwargs['outputdir']), kwargs['corpusprefix'], doc.id + ".salt")
                        })
                        # *metadata
                        )


    outputfile = os.path.join(kwargs['outputdir'], kwargs['corpusprefix'], doc.id + ".salt")
    xml = lxml.etree.tostring(saltdoc, xml_declaration=True, pretty_print=True, encoding='utf-8')
    with open(outputfile,'wb') as f:
        f.write(xml)
    print(f"Wrote {outputfile}",file=sys.stderr)
    return saltdoc, corpusdoc

def convert_tokens(doc, layers, **kwargs):
    """Convert FoLiA tokens (w) to salt nodes. This function also extracts the text layer and links to it."""
    audio_uri = "/Users/parkerglenn/Desktop/Brandeis/birch/folia2annis/010181.m4a"
    audio = True

    tokens = []
    textrelations = []
    phonrelations = []

    nodes_seqnr = 0

    map_id_to_nodenr = {}

    text = ""
    prevword = None

    phon = ""

    #will be initialised on first iteration:
    token_namespace = None
    layer = None

    if doc.declared(folia.TextContent):
        textnode = nodes_seqnr
        nodes_seqnr += 1
    if doc.declared(folia.PhonContent):
        if list(doc.select(folia.PhonContent)): #(not very efficient)
            phonnode = nodes_seqnr
            nodes_seqnr += 1

    # PARKER CHANGES
    # Create utterance node

    sstructure_ix = 0
    dom_rel_ix = 0 # Tracks id of current dominance relation ids

    """
    Tracks index of current node. Starts at 3, since birch salt documents begin with:
        - SAudioDS
        - STimeline
        - STextualDS
    """
    # node_ix = 2 if audio else 1
    node_ix = len(BASE_LAYERS)
    sstructure_data = OrderedDict() # Tracks id as keys, with node_ix and audio span info as values
    sspan_relations = OrderedDict() # Tracks word id as keys, with SSpan ids as values

    for word in doc.words():
        # print([i for i in word.json()['children'] if i['type'] == 'pos'])

        if word.parent.id not in sstructure_data:
            # Create Sstructure and link audio
            assert word.parent.json()['type'] == 'utt'
            tokens.append(create_sstructure(doc, sstructure_ix, **kwargs))
            sstructure_data[word.parent.id] = {}
            sstructure_data[word.parent.id]['node_ix'] = node_ix
            sstructure_data[word.parent.id]['audio_begin'] = time_tuple_to_string(word.parent.begintime)
            sstructure_data[word.parent.id]['audio_end'] = time_tuple_to_string(word.parent.endtime)
            sstructure_ix+=1
            node_ix += 1

        # Word is just the text token
        if not word.id:
            raise Exception("Only documents in which all words have IDs can be converted. Consider preprocessing with foliaid first.")
        if token_namespace is None:
            #only needs to be done once
            layer, token_namespace = init_layer(layers, word)
        textstart = len(text)
        try:
            text += word.text()
        except folia.NoSuchText:
            pass
        textend = len(text)

        phonstart = len(phon)
        try:
            phon += word.phon()
        except folia.NoSuchPhon:
            pass
        phonend = len(phon)

        word.nodes_seqnr = nodes_seqnr #associate it with the folia temporarily for a quick lookup later
        map_id_to_nodenr[word.id] = nodes_seqnr
        layer['nodes'].append(nodes_seqnr)
        # print("NODE SEQNR:")
        # print(nodes_seqnr)
        # print(node_ix)
        tokens.append(
            E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SToken",
                    },
                    *convert_identifier(word, **kwargs),
                    # *convert_type_information(word, **kwargs),
                    *convert_common_attributes(word, token_namespace, **kwargs),
                    *convert_inline_annotations(word, layers, **kwargs)
            )
        )
        # Node ids correspond to each items index in tokens list
        if word.parent.json()['type'] == 'utt':
            dom_rel_ix += 1
            tokens.append(create_dominance_relation(doc, dom_rel_ix, sstructure_data[word.parent.id]['node_ix'], node_ix, **kwargs))


        if text and textstart != textend:
            # print(word.text())
            # print(node_ix)
            # print(textnode)
            # print()
            textrelations.append(
                E.edges({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:STextualRelation",
                        "source": f"//@nodes.{node_ix}",
                        "target": f"//@nodes.{textnode}"
                        },
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                            "namespace": "salt",
                            "name": "id",
                            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sTextRel' + str(len(textrelations)+1)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SNAME",
                            "value": "T::sTextRel" + str(len(textrelations)+1)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SSTART",
                            "value": f"N::{textstart}"
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SEND",
                            "value": f"N::{textend}"
                        })
                )
            )

            if word.space or (prevword and word.parent != prevword.parent):
                text += " "

        if phon and phonstart != phonend:
            phonrelations.append(
                E.edges({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:STextRelation",
                        "source": f"//@nodes.{nodes_seqnr}",
                        "target": "//@nodes.{phonnode}"
                        },
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                            "namespace": "salt",
                            "name": "id",
                            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sPhonRel' + str(len(phonrelations)+1)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SNAME",
                            "value": "T::sPhonRel" + str(len(phonrelations)+1)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SSTART",
                            "value": f"N::{phonstart}"
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SEND",
                            "value": f"N::{phonend}"
                        })
                )
            )

        nodes_seqnr += 1
        prevword = word
        node_ix += 1

        """
        Birch specific: grab desc annotations in <pos></pos> and convert to SSpan.
        """
        word_ix = node_ix - 1
        pos_annotations = [i for i in word.json()['children'] if i['type'] == 'pos']
        for annotation in pos_annotations:
            description = [i for i in annotation['children'] if i['type'] == 'desc']
            if description:
                for feat in description[0]['value'].split(','):
                    feat = feat.strip()
                    if not feat.strip() in FEAT_DICT:
                        logger.yellow(f"ERROR: Unkown feature {feat}")
                    else:
                        sspan_id = str(sum([len(i) for i in sspan_relations.values()]))
                        tokens.append(
                            E.nodes({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpan",
                            },
                            E.labels({
                                 "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                "namespace": "salt",
                                "name": "id",
                                "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sSpan' + sspan_id
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                "namespace": "salt",
                                "name": "SNAME",
                                "value": "T::SSpan" + sspan_id
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                                    "namespace": "morph",
                                    "name": FEAT_DICT[feat],
                                    "value": "T::" + feat
                            }),
                        ))
                        # now we log to sspan_relations
                        if word_ix not in sspan_relations:
                            sspan_relations[word_ix] = [node_ix]
                        else:
                            sspan_relations[word_ix].append(node_ix)
                        node_ix += 1

    # Backtracking over accumulated sstructures outside of loop
    if audio:
        for ix, (sstructure_id, data) in enumerate(sstructure_data.items()):
            # print("ID:")
            # print(sstructure_id)
            textrelations.append(
                E.edges({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SAudioRelation",
                        "source": f"//@nodes.{data['node_ix']}",
                        "target": f"//@nodes.2"
                        },
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                            "namespace": "salt",
                            "name": "id",
                            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sAudioRel' + str(ix)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SNAME",
                            "value": "T::sAudioRel" + str(ix)
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SSTART",
                            "value": data['audio_begin']
                        }),
                        E.labels({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                            "namespace": "salt",
                            "name": "SEND",
                            "value": data['audio_end']
                        })
                )
            )


    # Backtracking over accumulated sspans and linking with SSpanningRelations
    sspan_rel_id=0
    for word_node_ix, sspan_node_ixs in sspan_relations.items():
        for sspan_node_ix in sspan_node_ixs:
            textrelations.append(
                 E.edges({
                     "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpanningRelation",
                         "source": f"//@nodes.{sspan_node_ix}",
                         "target": f"//@nodes.{word_node_ix}",
                     },
                 E.labels({
                     "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                    "namespace": "salt",
                    "name": "id",
                    "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#SSpanningRelation' + str(sspan_rel_id)
                }),
                E.labels({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                    "namespace": "salt",
                    "name": "SNAME",
                    "value": "T::SSpanningRelation" + str(sspan_rel_id)
                }),

            ))
            sspan_rel_id += 1


    nodes = []
    if text:
        if "sDocumentStructure:STextualDS" in BASE_LAYERS:
            nodes.append( E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:STextualDS",
                                },
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                    "namespace": "saltCommon",
                                    "name": "SDATA",
                                    "value": "T::" + text, #this can be huge!
                                }),
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                    "namespace": "salt",
                                    "name": "id",
                                    "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#TextContent'
                                }),
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                    "namespace": "salt",
                                    "name": "SNAME",
                                    "value": "T::TextContent"
                                }),
                        ))
    if audio:
        if "sDocumentStructure:STimeline" in BASE_LAYERS:
            nodes.append( E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:STimeline",
                                },
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                    "namespace": "salt",
                                    "name": "id",
                                    "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#TimelineContent'
                                }),
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                    "namespace": "salt",
                                    "name": "SNAME",
                                    "value": "T::TimelineContent"
                                }),
                        ))
        if "sDocumentStructure:SAudioDS" in BASE_LAYERS:
            nodes.append( E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SAudioDS",
                                },
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:sFeature",
                                    "namespace": "saltCommon",
                                    "name": "SAUDIO_REFERENCE",
                                    "value": "U::file:" + audio_uri
                                }),
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                    "namespace": "salt",
                                    "name": "id",
                                    "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#AudioContent'
                                }),
                                E.labels({
                                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                    "namespace": "salt",
                                    "name": "SNAME",
                                    "value": "T::AudioContent"
                                }),
                        ))

    if phon:
        nodes.append(E.nodes({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:STextualDS",
                            },
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                "namespace": "saltCommon",
                                "name": "SDATA",
                                "value": "T::" + phon, #this can be huge!
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                "namespace": "salt",
                                "name": "id",
                                "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#PhonContent'
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                "namespace": "salt",
                                "name": "SNAME",
                                "value": "T::PhonContent"
                            }),
                    ))
    nodes += tokens
    return (nodes, textrelations, text, phonrelations, phon, map_id_to_nodenr, nodes_seqnr, node_ix)

def init_layer(layers, element):
    """Initialises a salt layer (in a temporary structure) and computes the namespace for a certain annotation type"""
    if element.ANNOTATIONTYPE is None:
        raise Exception("Unable to init layer for element " + repr(element))
    namespace = "FoLiA::" + folia.annotationtype2str(element.ANNOTATIONTYPE).lower() + ("::" + element.set if element.set else "")
    if namespace not in layers:
        layers[namespace] = {
            "type": folia.annotationtype2str(element.ANNOTATIONTYPE).lower(),
            "set": element.set if element.set else "",
            "nodes": [],
            "edges": [],
        }
    return (layers[namespace], namespace)

def build_layers(layers, nodes, edges, node_ix):
    """Builds the final salt layers from the temporary structure and modifies nodes and edges accordingly (adding the @layer attribute)"""
    for i, (namespace, layer) in enumerate(layers.items()):
        for n in layer["nodes"]:
            # We don't want layer annotation on certain nodes
            if etree_to_dict(nodes[n])['@{http://www.w3.org/2001/XMLSchema-instance}type'] not in BASE_LAYERS:
                if n < len(nodes):
                    if "layers" in nodes[n].attrib:
                        nodes[n].attrib["layers"] += f" //@layers.{i}"
                    else:
                        nodes[n].attrib["layers"] = f"//@layers.{i}"
                if n < len(edges):
                    if "layers" in edges[n].attrib:
                        edges[n].attrib["layers"] += f" //@layers.{i}"
                    else:
                        edges[n].attrib["layers"] = f"//@layers.{i}"

        yield E.layers({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SLayer",
                    "nodes": " ".join(f"//@nodes.{n}" for n in range(2, node_ix+1)),
                    "edges": " ".join(f"//@edges.{n}" for n in layer["edges"])
                },
                E.labels({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                    "namespace": "salt",
                    "name": "id",
                    "value": "T::" + namespace
                }),
                E.labels({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                    "namespace": "salt",
                    "name": "SNAME",
                    "value": "T::" + layer['type'] + (" ("+ layer['set']+")" if layer['set'] else "")
                })
            )


def convert_inline_annotations(word, layers, **kwargs):
    """Convert FoLiA inline annotations to salt SAnnotation labels on a token"""
    for annotation in word.select(folia.AbstractInlineAnnotation):
        # if kwargs['saltnamespace'] and (annotation.set is None or annotation.set == word.doc.defaultset(annotation.ANNOTATIONTYPE)):
           #add a simplified annotation in the Salt namespace, this facilitates
           #conversion to other formats
       yield E.labels({
           "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                   "namespace": "salt",
                   "name": annotation.XMLTAG,
                   "value": "T::" + annotation.cls
               })
    """
    Birch specific: grab desc annotations in <pos></pos> and convert to SAnnotations.
    """
    # pos_annotations = [i for i in word.json()['children'] if i['type'] == 'pos']
    # for annotation in pos_annotations:
    #     description = [i for i in annotation['children'] if i['type'] == 'desc']
    #     if description:
    #         for feat in description[0]['value'].split(','):
    #             feat = feat.strip()
    #             if not feat.strip() in FEAT_DICT:
    #                 logger.yellow(f"ERROR: Unkown feature {feat}")
    #             else:
    #                 yield E.labels({
    #                     "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
    #                         "namespace": "salt",
    #                         "name": FEAT_DICT[feat],
    #                         "value": "T::" + feat
    #                 })


        # if not kwargs['saltonly']:
        #     layer, namespace = init_layer(layers, annotation)
        #     if word.nodes_seqnr is not None:
        #         layer['nodes'].append(word.nodes_seqnr)
        #
        #     for x in convert_common_attributes(annotation, namespace, **kwargs):
        #         yield x
        #
        #     for x in convert_features(annotation, namespace, **kwargs):
        #         yield x
        #
        #     for x in convert_higher_order(annotation, namespace, **kwargs):
        #         yield x

def convert_structure_annotations(doc, layers, map_id_to_nodenr, nodes_seqnr, **kwargs):
    """Convert FoLiA structure annotations (sentences, paragraphs, etc) to salt SSpan nodes and SSpaningRelation edges.
    In this conversion the structure annotations directly reference the underlying tokens, rather than other underlying structural elements like FoLiA does.
    """

    structure_nodes = []
    structure_spanningrelations = []
    #Create spans and text relations for all structure elements
    for structure in doc.select(folia.AbstractStructureElement):
        if not isinstance(structure, folia.Word): #word are already covered
            span_nodes = [ map_id_to_nodenr[w.id] for w in structure.words() ]
            if span_nodes:
                layer, namespace = init_layer(layers, structure)
                structure.nodes_seqnr = nodes_seqnr #associate it with the folia temporarily for a quick lookup later
                layer['nodes'].append(nodes_seqnr)
                structure_nodes.append(
                        E.nodes({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpan",
                            },
                            E.labels({ #we follow the example of the TCF->Salt converter here where it is used for sentences, a bit of weird entry, hopefully they knew what they were doing and this triggers some special behaviour for some of the converters? I just made it a bit more generic so it works for all structure types.
                                      "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                                "name": folia.annotationtype2str(structure.ANNOTATIONTYPE).lower(),
                                "value": "T::" + folia.annotationtype2str(structure.ANNOTATIONTYPE).lower()
                            }),
                            *convert_identifier(structure, **kwargs),
                            # *convert_type_information(structure, **kwargs),
                            *convert_common_attributes(structure, namespace, **kwargs),
                            *convert_features(structure, namespace, **kwargs),
                            *convert_higher_order(structure, namespace, **kwargs),
                        )
                )

                for nodenr in span_nodes:
                    structure_spanningrelations.append(
                         E.edges({
                             "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpanningRelation",
                                 "source": f"//@nodes.{nodes_seqnr}", #the structure
                                 "target": f"//@nodes.{nodenr}", #the token in the span
                             },
                             E.labels({
                                 "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                "namespace": "salt",
                                "name": "id",
                                "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sStructureSpanRel' + str(len(structure_spanningrelations)+1)
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                "namespace": "salt",
                                "name": "SNAME",
                                "value": "T::sStructureSpanRel" + str(len(structure_spanningrelations)+1)
                            }),
                        )
                    )
                nodes_seqnr += 1
    return (structure_nodes, structure_spanningrelations, nodes_seqnr)

def convert_span_annotations(doc, layers, map_id_to_nodenr, nodes_seqnr, **kwargs):
    """Convert FoLiA span annotations (sentences, paragraphs, etc) to salt SSpan nodes and SSpaningRelation edges.
    In this conversion the span annotations directly reference the underlying tokens, rather than other underlying structural elements like FoLiA does.
    """

    span_nodes = []
    span_spanningrelations = []
    #Create spans and text relations for all span elements:
    #only handles simple span elements that do not take span roles
    for span in doc.select(folia.AbstractSpanAnnotation):
        if not isinstance(span, (folia.AbstractSpanRole, folia.SyntacticUnit)) and  not any((isinstance(x, folia.AbstractSpanRole) for x in span.ACCEPTED_DATA)):
            span_token_nodes = [ map_id_to_nodenr[w.id] for w in span.wrefs() ]
            if span_token_nodes:
                layer, namespace = init_layer(layers, span)
                layer['nodes'].append(nodes_seqnr)

                span_nodes.append(
                        E.nodes({
                            "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpan",
                            },
                            *convert_identifier(span, **kwargs),
                            # *convert_type_information(span, **kwargs),
                            *convert_common_attributes(span, namespace, **kwargs),
                            *convert_features(span, namespace, **kwargs),
                            *convert_higher_order(span, namespace, **kwargs),
                        )
                )

                for nodenr in span_token_nodes:
                    span_spanningrelations.append(
                         E.edges({
                             "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SSpanningRelation",
                                 "source": f"//@nodes.{nodes_seqnr}", #the span
                                 "target": f"//@nodes.{nodenr}", #the token in the span
                             },
                             E.labels({
                                 "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                                "namespace": "salt",
                                "name": "id",
                                "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#sSpanRel' + str(len(span_spanningrelations)+1)
                            }),
                            E.labels({
                                "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                                "namespace": "salt",
                                "name": "SNAME",
                                "value": "T::sSpanRel" + str(len(span_spanningrelations)+1)
                            })
                        )
                    )
                nodes_seqnr += 1
    return (span_nodes, span_spanningrelations, nodes_seqnr)

def convert_syntax_annotation(doc, layers, map_id_to_nodenr, nodes_seqnr, **kwargs):
    syntax_nodes = []
    syntax_relations = []
    for syntaxlayer in doc.select(folia.SyntaxLayer):
        for su in syntaxlayer.select(folia.SyntacticUnit, recursive=False):
            nodes, relations, nodes_seqnr = convert_nested_span(su, layers, map_id_to_nodenr, nodes_seqnr, **kwargs)
            syntax_nodes += nodes
            syntax_relations += relations
    return (syntax_nodes, syntax_relations, nodes_seqnr)

def convert_nested_span(span, layers, map_id_to_nodenr, nodes_seqnr, **kwargs):
    nested_nodes = []
    nested_relations = []
    #process children first
    children_nodenr = []
    for child in span.select( (span.__class__, folia.Word), recursive=False):
        if child.__class__ is span.__class__:
            nodes, relations, nodes_seqnr = convert_nested_span(child, layers, map_id_to_nodenr, nodes_seqnr, **kwargs)
            nested_nodes += nodes
            nested_relations += relations
            children_nodenr.append(child.nodes_seqnr)
        else:
            children_nodenr.append(child.nodes_seqnr)

    layer, namespace = init_layer(layers, span)
    span.nodes_seqnr = nodes_seqnr #we will need this to get the node number back from the parent elements
    layer['nodes'].append(nodes_seqnr)

    nested_nodes.append(
            E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SStructure", #salt calls nested hierarchies 'structure', not to be confused with what FoLiA calls structure (document structure)
                },
                *convert_identifier(span, **kwargs),
                # *convert_type_information(span, **kwargs),
                *convert_common_attributes(span, namespace, **kwargs),
                *convert_features(span, namespace, **kwargs),
                *convert_higher_order(span, namespace, **kwargs),
            )
    )

    for nodenr in children_nodenr:
        nested_relations.append(
             E.edges({
                 "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SDominanceRelation",
                     "source": f"//@nodes.{nodes_seqnr}", #the span
                     "target": f"//@nodes.{nodenr}", #the subspan or token
                 },
                 E.labels({
                     "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                    "namespace": "salt",
                    "name": "id",
                    "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + span.doc.id + f"#sDomRel{nodes_seqnr}-{nodenr}"
                }),
                E.labels({
                    "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                    "namespace": "salt",
                    "name": "SNAME",
                    "value": f"T::sDomRel{nodes_seqnr}-{nodenr}"
                })
            )
        )
    nodes_seqnr += 1

    return (nested_nodes, nested_relations, nodes_seqnr)


def convert_identifier(annotation, **kwargs):
    if annotation.id:
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
            "namespace": "salt",
            "name": "id",
            "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + annotation.doc.id + '#' + annotation.id
        })
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
            "namespace": "salt",
            "name": "SNAME",
            "value": "T::" + annotation.id
        })

def create_sstructure(doc, sstructure_ix, **kwargs):
    """
    Takes in doc and cat, and creates sstructure node to connect others to.
    """
    return E.nodes({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SStructure"
        },
        E.labels({
        "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
        "namespace": "salt",
        "name": "id",
        "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#' + str(sstructure_ix)
        }),
        E.labels({
           "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                   "namespace": "salt",
                   "name": "cat",
                   "value": "T::" + 'UTT'
               }),
        E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
            "namespace": "saltCommon",
            "name": "SNAME",
            "value": "T::structure" + str(sstructure_ix), #this can be huge!
        })
    )

def create_dominance_relation(doc, dom_rel_ix, source_node_ix, target_node_ix, **kwargs):
    """
    Creates a dominance relation between a root_node_ix (pointing to a sstructure) and a connected_node_ix (probably pointing to an SToken).
    """
    return E.edges({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "sDocumentStructure:SDominanceRelation",
            "source": f"//@nodes.{source_node_ix}",
            "target": f"//@nodes.{target_node_ix}"
        },
        E.labels({
        "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
        "namespace": "salt",
        "name": "id",
        "value": "T::salt:/" + kwargs['corpusprefix'] + "/" + doc.id + '#' + str(dom_rel_ix)
        }),
        E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
            "namespace": "saltCommon",
            "name": "SNAME",
            "value": "T::structure" + str(dom_rel_ix), #this can be huge!
        })
    )




# def convert_type_information(annotation, **kwargs):
#     """Adds extra FoLiA type information as SMetaAnnotation labels (on salt nodes)"""
#     if annotation.XMLTAG:
#         yield E.labels({
#             "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SMetaAnnotation",
#             "namespace": "FoLiA",
#             "name": "elementtype",
#             "value": "T::" + annotation.XMLTAG
#         })
#     if annotation.ANNOTATIONTYPE:
#         yield E.labels({
#             "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SMetaAnnotation",
#             "namespace": "FoLiA",
#             "name": "annotationtype",
#             "value": "T::" + folia.annotationtype2str(annotation.ANNOTATIONTYPE).lower()
#         })

def convert_common_attributes(annotation,namespace, **kwargs):
    """Convert common FoLiA attributes as salt SMetaAnnotation labels (on salt nodes)"""

    if annotation.cls is not None:
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                    "namespace": namespace,
                    "name": "class",
                    "value": "T::" + annotation.cls
                })


def convert_features(annotation, namespace, **kwargs):
    """Convert FoLiA features to SAnnotation labels (on salt nodes)"""
    for feature in annotation.select(folia.Feature, recursive=False):
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SAnnotation",
                    "namespace": namespace,
                    "name": "feature/" + feature.subset,
                    "value": "T::" + feature.cls
                })

def convert_higher_order(annotation, namespace, **kwargs):
    """Convert certain FoLiA higher-order features to SAnnotation labels (on salt nodes)"""
    for seqnr, description in enumerate(annotation.select(folia.Description, recursive=False)):
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SMetaAnnotation",
                    "namespace": namespace,
                    "name": "description/" + str(seqnr + 1),
                    "value": "T::" + description.value
                })

    for seqnr, comment in enumerate(annotation.select(folia.Comment, recursive=False)):
        yield E.labels({
            "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SMetaAnnotation",
                    "namespace": namespace,
                    "name": "comment/" + str(seqnr + 1),
                    "value": "T::" + comment.value
                })


def convert_corpus(corpusdocs, **kwargs):
    edges = [ E.edges({
        "{http://www.w3.org/2001/XMLSchema-instance}type": "sCorpusStructure:SCorpusDocumentRelation",
                        "source": f"//@nodes.0",
                        "target": "//@nodes." + str(i+1),
                      },
                    E.labels({ # document ID
                              "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                        "namespace": "salt",
                        "name": "id",
                        "value": "T::salt:/corpDocRel" + str(i+1)
                    }),
                    E.labels({
                        "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                        "namespace": "salt",
                        "name": "SNAME",
                        "value": "T::corpDocRel" + str(i+1)
                    })
             )
             for i,_ in enumerate(corpusdocs) ]
    saltproject = getattr(E,"{saltCommon}SaltProject")( #sDocumentGraph
        {"{http://www.omg.org/XMI}version":"2.0"},
        E.sCorpusGraphs(
            E.labels({ # document ID
              "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                "namespace": "salt",
                "name": "SNAME",
                "value": "T::" +  kwargs['corpusprefix']
            }),
            E.nodes({
                "{http://www.w3.org/2001/XMLSchema-instance}type": "sCorpusStructure:SCorpus",
                    },
                    E.labels({ # document ID
                              "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SElementId",
                        "namespace": "salt",
                        "name": "id",
                        "value": "T::salt:/" + kwargs['corpusprefix']
                    }),
                    E.labels({
                        "{http://www.w3.org/2001/XMLSchema-instance}type": "saltCore:SFeature",
                        "namespace": "salt",
                        "name": "SNAME",
                        "value": "T::" + kwargs['corpusprefix']
                    }),
            ),
            *corpusdocs,
            *edges,
        )
    )
    outputfile = os.path.join(kwargs['outputdir'], "saltProject.salt")
    xml = lxml.etree.tostring(saltproject, xml_declaration=True, pretty_print=True, encoding='utf-8')
    with open(outputfile,'wb') as f:
        f.write(xml)
    print(f"Wrote project file {outputfile}",file=sys.stderr)
    return saltproject




def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','-V','--version',help="Show version information", action='version', version="FoLiA-tools v" + TOOLVERSION + ", using FoLiA v" + folia.FOLIAVERSION + " with library FoLiApy v" + folia.LIBVERSION, default=False)
    parser.add_argument('--recurse','-r',help="recurse into subdirectories", action='store_true', required=False)
    parser.add_argument('--extension','-e',type=str, help="extension", action='store', default="xml",required=False)
    parser.add_argument('--corpusprefix','-p', type=str, help="Corpus prefix for salt", action='store', default="foliacorpus",required=False)
    parser.add_argument('--saltnamespace','-s',help="Add simplified annotations in the salt namespace", action='store_true', required=False, default=False)
    parser.add_argument('--outputdir','-o',type=str, help="Output directory", action='store', default=".", required=False)
    parser.add_argument('--saltonly','-S',help="Skip complex annotations not in the salt namespace (use with --saltnamespace). This will ignore a lot of the information FoLiA provides!", action='store_true', required=False, default=True)
    parser.add_argument('files', nargs='*', help='Files (and/or directories) to convert')
    args = parser.parse_args()

    os.makedirs(os.path.join(os.path.realpath(args.outputdir), args.corpusprefix), exist_ok=True)
    kwargs = args.__dict__
    if args.files:
        corpusdocs = []
        skipnext = False
        for file in args.files:
            if os.path.isdir(file):
                corpusdocs += processdir(file, **kwargs)
            elif os.path.isfile(file):
                saltdoc, corpusdoc = convert(file, **kwargs)
                corpusdocs.append(corpusdoc)
            else:
                print("ERROR: File or directory not found: " + file,file=sys.stderr)
                sys.exit(3)
        if not corpusdocs:
            sys.exit(1)
        convert_corpus(corpusdocs, **kwargs)
    else:
        print("ERROR: No files specified. Add --help for usage details.",file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()



def testing():
    import pickle
    with open('/Users/parkerglenn/Desktop/Brandeis/birch/foliatools/foliatools/example_args.pkl', 'rb') as f:
        args = pickle.load(f)
    kwargs = args.__dict__
    file = '/Users/parkerglenn/Desktop/Brandeis/birch/folia-testing/folia_test_document/pos.2.0.0.folia.xml'

    doc = convert(file, **args.__dict__)

    dir(doc)
    for word in doc.words():
        print(word)

    s = word.parent
    dir(s)
    word.parent.json()['type']
    s.json()
    s.id
    s.xmlstring()

    s.text()
    dir(word)


    word, token_namespace = token_namespaces[0]
    dir(word)
    word.elements()
    for annotation in word.select(folia.AbstractInlineAnnotation):
        print(annotation.cls)



    for word in [i[0] for i in token_namespaces]:
        print(word.text())
"""
TODO:
    - Bump up STextualRelation IDs by 1, since we're adding in SStructures that increment the node count
"""
