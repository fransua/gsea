#!/usr/bin/python
"""
17 Aug 2010


"""

__author__     = "Francois-Jose Serra"
__email__      = "francois@barrabin.org"
__licence__    = "GPLv3"
__version__    = "0.1b"
__title__      = "gene set tool kit v%s" % __version__
__references__ = '''
        Al-Shahrour, F., Diaz-Uriarte, R., & Dopazo, J. 2004.
             FatiGO: a web tool for finding significant associations of Gene Ontology terms with groups of genes.
             Bioinformatics (Oxford, England) 20: 578-80.
             Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/14990455
        Al-Shahrour, F., Diaz-Uriarte, R., & Dopazo, J. 2005.
             Discovering molecular functions significantly related to phenotypes by combining gene expression data and biological information.
             Bioinformatics (Oxford, England) 21: 2988-93.
             Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/15840702
'''


# easy_install fisher
try:
    from fisher import pvalue
except ImportError:
    # or compile it somewhere...
    # http://pypi.python.org/pypi/fisher
    from extra_stats.fisher import pvalue
# in my extra_stats package
from stats.fdr import bh_qvalues
from numpy import log
from sys import stderr
from optparse import OptionParser
from bisect import bisect_left
from warnings import warn

class GSEA:
    '''
    Fatiscan with upper case, it is an object.
    USAGE:

      from gsea import GSEA
      infile = "data_test/Homo_sapiens_1.val"
      annot  = "data_test/biol_proc_2-8.annot"
      gene_set = GSEA (infile, annot, use_order=False)
      list1 = gene_set.genes[10:250]
      list2 = gene_set.genes[900:1500]
      #gene_set.run_fatigo (list1)
      gene_set.run_fatiscan()
      
      gene_set.summarize("GO:0000278")
      
    '''

    def __init__ (self, gene_vals, annot, algo='fatiscan', use_order=True):
        '''
        init function, what is done when object is called.
        '''
        # get gene list and corresponding values
        self.use_order = use_order

        if algo == 'fatiscan':
            self.genes, self.values, self.order = self._parse_values (gene_vals)
            self.annot = self._parse_annot(annot)
            # sort genes in annot by their values...
            # useful to know which gene in list1 or list2
            self.annot   = self._order_genes_in_annot()
        elif algo == 'fatigo':
            self.annot = self._parse_annot(annot)
        self.gsea    = {}
        self._thresh = []

    def _parse_values (self, gene_vals):
        '''
        parse in file in format:
        geneID_1 <tab> value1
        geneID_2 <tab> value2
        ...
        genes should be ordered by value
        returns genes, values and order of values
        '''
        if type (gene_vals) is list :
            genes, values = zip (* sorted (gene_vals))
        else:
            genes, values = zip (*sorted ((i.strip().split('\t')
                                           for i in open(gene_vals)),
                                          key=lambda x: float(x[1])))
        values = map (float, values)
        if self.use_order:
            order = map (values.index, values)
        else:
            order = values[:]
        return genes, dict (zip (genes, values)), order

    def _parse_annot (self, annot):
        '''
        parse annotation file in format:
        annotationA <tab> geneID_1
        annotationA <tab> geneID_2
        ...
        speed notes: * iterator on for for
                     * dico
        '''
        dico = {}
        if type (annot) is dict :
            dico = annot
            return dico
        if hasattr (self, 'values'):
            # in case we already have a list of genes we can reduce annotation
            # to our genes
            values = self.values
            for line in open (annot):
                gene, annot = line.strip().split('\t')
                if values.has_key (gene):
                    dico.setdefault (annot, []).append (gene)
        else:
            for line in open (annot):
                gene, annot = line.strip().split('\t')
                dico.setdefault (annot, []).append (gene)
        return dico

    def _order_genes_in_annot(self):
        '''
        order genes in annot dict by their values
        '''
        dico = {}
        for annot in self.annot.iterkeys():
            dico[annot] = sorted (self.annot[annot], \
                                  key=lambda x: self.values[x])
        return dico

    def _adjust_pvals (self):
        '''
        adjusts pvalues of a list of fatigos using fdr
        '''
        pvalues = []
        for nam in sorted (self.gsea):
            pvalues.append (self.gsea[nam]['pv'])
        qvalues = iter (bh_qvalues (pvalues))
        for nam in sorted (self.gsea):
            self.gsea[nam]['apv'] = qvalues.next()

    def _enrichment (self, genes1, genes2, part=''):
        '''
        computes an enrichment test between two sets of lists of genes
        '''
        len_genes1 = len (genes1)
        len_genes2 = len (genes2)
        dico = self.gsea
        part = '|' + part
        # start fishers
        for ann, annot_genes in self.annot.iteritems():
            p1 = len (set (annot_genes) & genes1)
            p2 = len (set (annot_genes) & genes2)
            n1 = len_genes1 - p1
            n2 = len_genes2 - p2
            dico [ann + part] = {'p1' : p1, 'n1': n1,
                                 'p2' : p2, 'n2': n2,
                                 'pv' : pvalue (p1, n1, p2, n2).two_tail, # 3/4 of time spent here
                                 'odd': _get_odd_ratio (p1, p2, n1, n2)}

    def run_fatigo (self, genes1, genes2=None):
        '''
        run fatigo analysis

        ref: Al-Shahrour, F., Diaz-Uriarte, R., & Dopazo, J. 2004.
             FatiGO: a web tool for finding significant associations of Gene Ontology terms with groups of genes.
             Bioinformatics (Oxford, England) 20: 578-80.
             Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/14990455
        '''
        genes1 = set (genes1)
        if not genes2:
            genes2 = set (self.genes) - set (genes1)
        else:
            genes2 = set (genes2)
        self._enrichment (genes1, genes2)
        self._adjust_pvals()
        def summarize(annot, with_part=False):
            '''
            returns result for most significant partition for given annot
            '''
            result = []
            result.append (self.gsea [annot + '|'])
            part, col = min (enumerate (result),
                             key=lambda (y,x): (x['apv'],y))
            if with_part:
                return part, col
            else:
                return col
        self.__dict__['summarize'] = summarize

    def run_fatiscan (self, partitions=30, verb=False):
        '''
        run gsea needs python fisher, and fdr from extra stats
        speed notes: * making annot genes and order local does not
                       significantly speed up process.
                     * putting external method with bissect part neither
                     * no more ideas...

        ref: Al-Shahrour, F., Diaz-Uriarte, R., & Dopazo, J. 2005.
             Discovering molecular functions significantly related to phenotypes by combining gene expression data and biological information.
             Bioinformatics (Oxford, England) 21: 2988-93.
             Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/15840702
        '''
        order_index = self.order.index
        self.gsea = {}
        rank = float (max (self.order) - min (self.order))/partitions
        # define cutoff value for each partition
        self._below_thresh = [bisect_left (self.order, rank * (part + 1)) \
                              for part in xrange(partitions)]
        self._thresh = [rank * (part + 1) for part in xrange(partitions)]
        all_genes = set (self.genes)
        for part in xrange(partitions):
            genes1 = set (self.genes [:self._below_thresh [part]])
            if verb:
                print 'partition %2s, thresh value= %-8s number of genes: %d'\
                      % (part, self._thresh [part], len (genes1))
            self._enrichment (genes1, all_genes - genes1, str (part))
        self._adjust_pvals ()
        def summarize(annot, with_part=False):
            '''
            returns result for most significant partition for given annot
            '''
            result = []
            for part in xrange (partitions):
                result.append (self.gsea [annot + '|' + str (part)])
            part, col = min (enumerate (result),
                             key=lambda (y,x): (x['apv'],y))
            if with_part:
                return part, col
            else:
                return col
        self.__dict__['summarize'] = summarize

    def write_gsea (self, outfile, max_apv=1, all_parts=False):
        '''
        write to file, or pickle
        '''
        def _get_string(dico, annot):
            '''
            get string from gsea_dic current value
            '''
            string = []
            string.append (dico ['p1' ])
            string.append (dico ['n1' ])
            string.append (dico ['p2' ])
            string.append (dico ['n2' ])
            string.append (dico ['odd'])
            string.append (dico ['pv' ])
            string.append (dico ['apv'])
            string.append (','.join (list (annot[:dico['p1']])))
            string.append (','.join (list (annot[dico['p1']:])))
            return map (str, string)
        if self.gsea == []:
            raise Exception ('ERROR: you do not have run GSEA yet...')
        cols = ['#term', 'part', 'list1_positives', 'list1_negatives',
                'list2_positives', 'list2_negatives', 'odds_ratio_log',
                'pvalue', 'adj_pvalue', 'list1_positive_ids (ordered by value)',
                'list2_positive_ids (ordered by value)']
        out = open (outfile, 'w')
        out.write('\t'.join(cols)+'\n')
        if all_parts:
            for annot in self.annot:
                for part in xrange (30):
                    nam = annot + '|' + str (part)
                    if self.gsea [nam]['apv'] > max_apv:
                        continue
                    string = _get_string (self.gsea[nam],
                                          self.annot[annot])
                    out.write ('\t'.join ([annot] + [str(part)] +string) + '\n')
        else:
            for annot in self.annot:
                part, col = self.summarize (annot, True)
                string = _get_string (col, self.annot [annot])
                out.write ('\t'.join ([annot] + [str(part)] + string) + '\n')
        out.close()

def main ():
    '''
    for direct command line call
    '''
    opts = get_options()
    gene_set = GSEA(opts.infile, opts.annot, algo=opts.algo,
                    use_order=opts.use_order)
    if opts.algo == 'fatiscan':
        gene_set.run_fatiscan (partitions=opts.partitions, verb=opts.verb)
    elif opts.algo == 'fatigo':
        gene_set.run_fatigo (map (lambda x:x.strip(),
                                  open (opts.list1).readlines()),
                             map (lambda x:x.strip(),
                                  open (opts.list2).readlines()))
    if opts.pickle:
        from cPickle import dump
        dump (open (outfile, 'w'), self)
    else:
        gene_set.write_gsea(opts.outfile, all_parts=opts.all_parts, \
                            max_apv=float(opts.max_apv))


def _get_odd_ratio (p1, p2, n1, n2):
    '''
    computes fisher odd ratio
    '''
    try:
        odd = log ((float (p1)/n1) / (float (p2)/n2))
    except ZeroDivisionError:
        if p2 == 0:
            odd = float ('inf')
        elif n1 == 0:
            warn ("WARNING: Empty partition(s)")
            odd = float ('-inf')
        elif n2 == 0:
            warn ("WARNING: Empty partition(s)")
            odd = float ('inf')
    return odd


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Gene set enrichment analysis                                                                                                
.                                                                           
********************************************                                      
"""
        )
    parser.add_option('-i', dest='infile', metavar="PATH",
                      help='''path to input file with a ranked list, in format:                             
                      geneID_1 <tab> value1                                                        
                      geneID_2 <tab> value2                                                         
                      ...
                      ''')
    parser.add_option('-x', dest='list1', metavar="PATH",
                      help='''path to first input file with a list of genes, in format:                             
                      geneID_1                                                        
                      geneID_2                                                         
                      ...
                      ''')
    parser.add_option('-y', dest='list2', metavar="PATH",
                      help='''path to second input file with a list of genes, in format:                             
                      geneID_1                                                        
                      geneID_2                                                         
                      ...
                      ''')
    parser.add_option('-a', dest='annot', metavar="PATH", \
                      help='''path to annotation file in format:                                           
                      annotationA <tab> geneID_1                                                       
                      annotationA <tab> geneID_2                                                   
                      ...
                      ''')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file tab separated file')
    parser.add_option('-R', '--use_rank', action='store_true', \
                      dest='use_order', default=False, \
                      help=\
                      '''[%default] Use rank of genes in stead of provided value. This option will smooth wired distributions and equalize partition sizes.''')
    parser.add_option('-p', metavar="INT", dest='partitions', default=30, \
                      help='''[%default] Number of partitions.''')
    parser.add_option('-F', dest='algo', default='fatiscan', \
                      help='''[%default] default algorithm.                              
                        * fatiscan needs -i                                         
                        * fatigo needs -x and -y ''')
    parser.add_option('--max_apv', metavar="FLOAT", dest='max_apv', default=1, \
                      help='''[%default] Only write to outfile results with adjusted pvalue higher than specified value.''')
    parser.add_option('--long', dest='all_parts', action='store_true', default=False, \
                      help='''[%default] Write results for all partitions.''')
    parser.add_option('--pickle', action='store_true', \
                      dest='pickle', default=False, \
                      help='[%default] Store results in python dict (cPickle) format.')
    parser.add_option('--verbose', action='store_true', \
                      dest='verb', default=False, \
                      help=\
                      '[%default] Talk a bit... ')
    opts = parser.parse_args()[0]
    if (not opts.infile and (not opts.list1 or not opts.list2)) \
           or not opts.annot or not opts.outfile:
        exit(parser.print_help())
    return opts




if __name__ == "__main__":
    exit(main())
