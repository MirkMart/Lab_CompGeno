import treeswift
import argparse

def unroot(tree):
    """
    Unroots treeswift tree. Adapted from treeswift 'deroot' function.
    This one doesn't contract (A,B); to A;

    Parameters
    ----------
    tree: treeswift tree

    Returns unrooted treeswift tree
    """
    if tree.root == None:
        return tree
    if tree.root.num_children() == 2:
        [left, right] = tree.root.child_nodes()
        if not right.is_leaf():
            right.contract()
        elif not left.is_leaf():
            left.contract()
    tree.is_rooted = False
    return tree


def reroot_on_edge(tree, node):
    if not node.is_root(): 
        if not hasattr(node, 'edge_length') or node.edge_length is None or node.edge_length == 0:
            node.edge_length = 1
        tree.reroot(node, length=node.edge_length / 2)


def remove_in_paralogs(tree, gene_to_species=lambda x:x):
    """
    Removes in-paralogs from unrooted tree.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns number of in-paralogs removed
    """
    # root tree if not rooted
    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])
        #tree.reroot(tree.root)

    num_paralogs = 0
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([gene_to_species(node.get_label())])
        else:
            node.s = set([])
            for child in node.child_nodes():
                node.s = node.s.union(child.s)
            # collapse if paralogs
            if len(node.s) == 1:
                for child in node.child_nodes()[1:]:
                    node.remove_child(child)
                    num_paralogs += 1

    # check over top of root
    if any(len(child.s) == 1 for child in tree.root.child_nodes()):
        for node in tree.traverse_preorder():
            if not node.is_root():
                parent = node.get_parent()
                node.up = set([]) if parent.is_root() else parent.up
                for sibl in parent.child_nodes():
                    if sibl != node:
                        node.up.union(sibl.s)
                if len(node.up) == 1:
                    for sibl in parent.child_nodes():
                        if sibl != node:
                            parent.remove_child(sibl)
                            num_paralogs += 1
    tree.suppress_unifurcations()
    return num_paralogs


def get_min_root(tree, gene_to_species=lambda x:x, verbose=False):
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns vertex corresponding to best edge to root tree on
    """

    def score(total_set, set1, set2):
        if not len(set1.intersection(set2)) == 0:
                if total_set == set1 or total_set == set2:
                    if set1 == set2:
                        return 1
                    else:
                        return 2
                else:
                    return 3 
        return 0

    # check if tree is single leaf
    if tree.root.num_children() == 0:
        tree.root.s = set([tree.root.get_label()])
        return tree.root, 0, [] 

    # root tree if not rooted
    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])
    tree.resolve_polytomies()

    # Get down scores pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([gene_to_species(node.get_label())])
            node.d_score = 0
        else:
            if node.num_children() != 2:
                raise Exception("Vertex has more than 2 children")

            [left, right] = node.child_nodes()
            node.down = left.down.union(right.down)
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down)

    min_score, best_root, ties = float("inf"), None, []

    # Get scores above edge pass
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            root.skip = True
        else:
            node.skip = False

    # Get 'up' set for children of root
    [left, right] = root.child_nodes()
    left.up = right.down
    left.u_score = right.d_score
    right.up = left.down
    right.u_score = left.d_score
    left.skip = True
    right.skip = True

    min_score = left.u_score + left.d_score + score(left.up.union(left.down), left.up, left.down)
    # we don't want to root at a leaf
    if not left.is_leaf():
        best_root = left 
    elif not right.is_leaf():
        best_root = right
    # if both are leaves (i.e. two leaf tree), we want to keep the same rooting
    else:
        best_root = root
    ties = [best_root]

    for node in tree.traverse_preorder(leaves=False):
        if not node.skip:
            
            parent = node.get_parent()
            if parent.child_nodes()[0] != node:
                other = parent.child_nodes()[0]
            else: 
                other = parent.child_nodes()[1]

            node.up = parent.up.union(other.down)
            node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down)

            total_score = node.u_score + node.d_score + score(node.up.union(node.down), node.up, node.down)

            if total_score == min_score:
                ties.append(node)
                
            if total_score < min_score:
                num_ties = 0
                min_score = total_score
                best_root = node
                ties = [node]

    if verbose:            
        print('Best root had score', min_score, 'there were', len(ties), 'ties.')
        
    return best_root, min_score, ties


def tag(tree, gene_to_species=lambda x:x):
    """
    Tags tree according to its current rooting.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label
    """
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([gene_to_species(node.get_label())])
            node.n_dup = 0
        else:
            [left, right] = node.child_nodes()

            node.s = left.s.union(right.s)
            node.n_dup = left.n_dup + right.n_dup
            if len(left.s.intersection(right.s)) == 0:
                node.tag = 'S'
            else: 
                node.tag = 'D'
                node.n_dup += 1
    tree.n_dup = tree.root.n_dup


def decompose(tree, single_tree=False):
    """
    Decomposes a tagged tree, by separating clades at duplication vetices

    NOTE: must be run after 'tag()'

    Parameters
    ----------
    tree: tagged treeswift tree
    single_tree: return only the single large

    Returns result of the decomposition as a list of trees
    """
    out = []
    for node in tree.traverse_postorder(leaves=False):
        if node.tag == 'D':
            # trim off smallest subtree (i.e. subtree with least species)
            [left, right] = node.child_nodes()
            delete = left if len(left.s) < len(right.s) else right
            if not single_tree:
                out.append(tree.extract_subtree(delete))
                out[-1].suppress_unifurcations()
            node.remove_child(delete)
    tree.suppress_unifurcations() # all the duplication nodes will be unifurcations
    out.append(tree)
    return out


def relabel(tree, gene_to_species=lambda x:x):
    for l in tree.traverse_postorder(internal=False):
        l.set_label(gene_to_species(l.get_label()))
    return tree


def main(args):
    
    if args.delimiter is not None:
        gene_to_species = lambda x : args.delimiter.join(x.split(args.delimiter)[:args.nth_delimiter])
    else:
        gene_to_species = lambda x : x

    if args.output is None:
        split = args.input.rsplit('.', 1)
        output = split[0] + '-decomp.' + split[1]
    else:
        output = args.output

    # delete existing outgroup file (so you don't append to it)
    outgroup_file_name = args.input.rsplit('.', 1)[0] + '_outgroups.txt'
    if args.outgroups:
        open(outgroup_file_name, 'w').close()

    with open(args.input, 'r') as fi, open(output, 'w') as fo:
        for i, line in enumerate(fi, 1):
            tree = treeswift.read_tree_newick(line)
            if type(tree) == list:
                assert len(tree) == 0, "Could not interpret {} on line {:d} as a tree".format(line, i)
                continue

            if args.remove_in_paralogs:
                num_paralogs = remove_in_paralogs(tree, gene_to_species)

            root, score, ties = get_min_root(tree, gene_to_species)
            reroot_on_edge(tree, root)
            tag(tree, gene_to_species)

            if args.verbose:
                print('Tree ', i, ': Tree has ', len(tree.root.s), ' species.', sep='')
                if args.remove_in_paralogs:
                    print(num_paralogs, 'in-paralogs removed prior to rooting/scoring.')  
                if len(tree.root.s) < 2:
                    print('Uninformative')
                elif tree.n_dup == 0:
                    print('Single-Copy')                 
                else:
                    outgroup = min((len(child.s), child.s) for child in tree.root.child_nodes())                    
                    print('Best root had score ', score, ' with ', tree.n_dup, ' non-terminal' if args.remove_in_paralogs else '',
                        ' duplications; there were ', len(ties), ' ties.\nOutgroup: {',','.join(outgroup[1]),'}', sep='')

            # Choose modes
            if args.no_decomp:
                out = [tree]
            else:
                out = list(filter(lambda x:x.num_nodes(internal=False) >= args.minimum, decompose(tree, args.single_tree)))

            # Output trees
            for t in out:
                if not args.no_decomp: unroot(t)
                if not args.keep_labels: relabel(t, gene_to_species)
                t.suppress_unifurcations()
                fo.write(t.newick() + '\n')
            
            if args.verbose:
                print('Decomposition strategy outputted', len(out), 'tree(s) with minimum size', args.minimum, '.\n')

            # output outgroups
            if args.outgroups:
                og_tree = treeswift.read_tree_newick(line)
                root, score, ties = get_min_root(og_tree, args.delimiter)
                reroot_on_edge(og_tree, root)
                tag(og_tree, args.delimiter)
                if len(og_tree.root.s) >= 2 and og_tree.n_dup >= 1:
                    with open(outgroup_file_name, 'a') as outgfile:
                        outgfile.write('Tree ' + str(i) + ':\n')
                        for t in ties:
                            reroot_on_edge(og_tree, t)
                            tag(og_tree, args.delimiter)
                            outgroup = min((len(child.s), child.s) for child in og_tree.root.child_nodes())
                            outgfile.write('{' + ','.join(outgroup[1]) + '}\n')
            

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='====================== DISCO v1.4.0 ======================')

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree list file")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")
    parser.add_argument('-n', '--nth-delimiter', type=int, # Default is 1 -- set below
                        help="Split on nth delimiter (only works with -d)")
    parser.add_argument('-m', "--minimum", type=int, 
                        help="Minimum tree size outputted", default=4)
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enables verbose output")
    parser.add_argument("--keep-labels", action='store_true', 
                        help="Keep original leaf labels instead of relabeling them with their species labels (only relevant with delimiter)")
    parser.add_argument('--single_tree', action='store_true',
                        help="Only output single large tree")
    parser.add_argument('--no-decomp', action='store_true', 
                        help="Outputs rooted trees without decomposition")
    parser.add_argument("--outgroups", action='store_true',
                        help="Output outgroups to file (including ties)")
    parser.add_argument("--remove_in_paralogs", action='store_true',
                        help="Remove in-paralogs before rooting/scoring tree.")

    args = parser.parse_args()
    if args.delimiter is None:
        if args.nth_delimiter is not None:
            parser.error("Cannot set -n without a delimiter")
        if args.keep_labels:
            parser.error("Cannot use --keep-labels without a delimiter")
    elif args.nth_delimiter is None:
        args.nth_delimiter = 1
    if args.single_tree and args.no_decomp:
        parser.error("Cannot combine --single_tree and --no-decomp")
    if not args.verbose and args.remove_in_paralogs:
        print("--remove_in_paralogs is meaningless without --verbose, as it does not change the optimal rooting. " + 
            "It may also slow the program.")

    main(args)
