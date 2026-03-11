def remove_tips(G):

    tips_removed = 0

    for node in list(G.nodes()):

        if G.out_degree(node) == 0 and G.in_degree(node) == 1:

            pred = list(G.predecessors(node))[0]

            if G.out_degree(pred) > 1:
                G.remove_node(node)
                tips_removed += 1

    return tips_removed


def remove_bubbles(G):

    bubbles_removed = 0

    for node in list(G.nodes()):

        if G.out_degree(node) > 1:

            succ = list(G.successors(node))

            if len(succ) == 2:

                s1, s2 = succ

                w1 = G[node][s1]["weight"]
                w2 = G[node][s2]["weight"]

                if w1 > w2:
                    G.remove_edge(node, s2)
                else:
                    G.remove_edge(node, s1)

                bubbles_removed += 1

    return bubbles_removed
