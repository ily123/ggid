"""Module to keep all of the app's text-heavy elements."""
# why not keep these in a csv file or database?


auc_explanation = {
    "random": """
        AUC ~= 0.50: information does not diffuse between the
        given kinases.
        This is because the kinases are not closely connected in
        the network.
        Either the network is not modeling their
        relationship properly, or these kinases are not
        very similar (i.e. do not share low-level GO terms).
        You should diffuse each kinase on its own
        to find its functional neighbors.
        """,
    "moderate": """
        0.60 < AUC < 0.80: there is some signal! The kinases in
        the set diffuse some information to each other,
        meaning they are somewhat connected in the network.

        Alternatively, there could
        be two independent kinase clusters of in your set.
        Kinases within the same cluster are predictive of each
        other, but not of the kinases in the other cluster.
        Each cluster may be responsible for some unique
        sub-function relevant to the overall process.
        Identify the clusters by looking at the graph visual
        above, and diffuse each cluster on its own to find
        its functional neighbors.
        """,
    "high": """
        AUC > 0.80: kinases are strongly predictive of each
        other. They are well connected in the network, and thier
        neighbors may be involved in the same type of process.
        For the list of their closest neigbors see table below.
        """,
}


# Lines run over the char limit -- just leave them
tabs = {
    "theory_tab": """
        ### Graphical Abstract
        ---
        ![graphical_abstract](assets/graphical_abstract.png)

        ### The problem: Using networks to annotate proteins
        ---
        Proteins don't work in isolation. They work together and form networks.

        In practice, this means that once you've identified a protein associated with some critical
        process, you can query that protein's network to find
        additional proteins that are involved in the same process.
        In case of disease, this allows you to expand the candidate list of theraupetic targets, and arrive at
        a cure faster.

        Construction of protein networks, and learning from them, is one of the
        subfields in computational biology. The industry standard for collecting and
        visualizing these data is the [STRING](https://string-db.org/cgi/about) database. In
        STRING
        you can find protein networks for many species, and for a variety of underlying data (physical
        interactions, high-throughput co-expression experiments, implied functional association, text mining, etc).

        In this dashboard, I present a new approach for creating protein networks that is
        complementary to STRING.
        I use data from the [Gene Ontology](https://geneontology.org/) database to
        create GO similarity networks. I then apply information diffusion to the network
        to find clusters of connected proteins. Below I explain the methodology.

        For a complete example for using the dashboard, see the EXAMPLE tab.

        ### Gene Ontology (GO) similarity network
        ---

        [Gene Ontology](http://geneontology.org/) database is to proteins what the [IMDB](https://www.imdb.com)
        is to movies. It's a database that catalogues proteins, and assigns them
        descriptive labels. Whereas IMDB may label movies as
        "comedy", "action", etc, the Gene Ontology annotates proteins
        with labels such as "protein kinase", "mitochondirial process", "cell wall",
        and so on. These labels are  known as "GO terms" and they describe a protein's
        job within the cell. Because proteins that do
        similar jobs have similar GO terms, GO term similarity can be readily
        converted into a network/graph (which we we can then learn from).

        To build a protein network from GO terms, we follow these steps:

        1. Given a set of proteins, we measure GO term similariy between each possible pair
         of proteins in the set. To measure GO term similarity between two proteins, we
         do the following:
            * We define GO term similarity of two proteins as the average similarity
            of their GO terms.
            * Similarity between two GO terms ```q``` and ```w``` is the specificity of
            their most informative common ancestor (```MICA```):

            ![equation](https://latex.codecogs.com/gif.latex?similarity(q,w)&space;=&space;specificity(MICA))

            * What is an ancestor GO term?

            GO terms are not just flat labels. GO terms
            are related to each other and are organized in a tree, where each term
            has ancestors and children. Ancestor terms are more general and high-level, while
            children terms are more specific and low-level. For example, ```Signaling``` is a
            general high-level term, and one of its children is the more specific ```Kinase signaling```.
            The child-ancestor relationships branch out with ancestors having multiple
            children.

             Thus, given two terms, we can trace their ancestry and find common ancestor
             terms.

            * Once we find the shared ancestor terms, how do we decide which is the
            most specific? That's easy. The specificity is simply the rarity of
            the GO term, eg how frequently the term appears in the Ontology:

            ![equation](https://latex.codecogs.com/gif.latex?specificity(a)&space;=&space;-log\\frac{n}{N})

            In the equation above, ```n``` is the number of times that term ```a``` has ever
            been annotated to a protein, while ```N``` is the total number of all annotations
            made. (Note that when we assign a GO term to a protein, we refer to it as an
            annotation, i.e. we annotate a protein with a GO term.) So ```n/N``` is simply the frequency of term in the annotation
            corpus. To make things more intuitive, we convert frequency into a positive
            number by taking a negative log of it (smaller frequency corresponds to higher negative log number).

            In practise, the most informative common ancestor is usually the lowest common ancestor in the tree.
            This is because terms lower in the tree are used less frequently, and are more specific.

            * To arrive at a final similarity score for two proteins, we measure the
            similarity of their GO terms (for every possible GO term pair), and take the
            Best Match Average (BMA) of the individual GO term similarity scores. (In BMA, rather
            than do a mean over all go term pairs, we instead select the highest score for each term
            and average those.) The BMA GO term similarity score is the similarity score for the two proteins.


        2. Once we have computed similarity for every protein pair in the set,
        we have a similarity matrix. In this matrix, higher scores correspond to pairs
        of highly-similar proteins, while low scores correspond to proteins that are
        dissimilar. To convert the matrix into a graph, we go trough each
            protein (ie column/row in the matrix) and set the top 5 of its most similar
            partner proteins as edges, and discrad all other pairs. This produced a protein
            graph where each protein is connected to 5 proteins that it is most similar to.

        3. To learn from this network, we then use information diffusion
        as described below.

        ### Graph Information Diffusion
        ---

        To extract information from the network, we can simply examine it by hand. Given a protein, what
        are its closest neighbors? This works, but in practise gets out of hand quickly
        when you want to examine more than a single protein at the same time, or when the number of edges
        is too high. So it's best to use an algorithmic approach, such as Information
        Diffusion.

        Intuitively, Information Diffusion is exactly what it sounds like. We label
        some nodes in the network with information, then let information diffuse to other
        nodes using edges as conduits. At the end of the diffusion experiment, nodes that retain the
        most information are the ones that are most
        closely connected to the original labeled nodes.

        Using this method, we can discover clusters of proteins that do the same job.
        For example, label all proteins associated with ```disease X```, diffuse
        that label, and measure amount of information in each node post diffusion.
        Rank the originally unlabled nodes by the amount of information they retained.
        Proteins at the top of the ranked list are likely to be involved in
        ```disease X```.

        For formal mathematic definition of diffusion see the methods section
        of [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3001439/).

        ---
        However, here is a quick summary from someone who isnt that great at linear
        algebra:

        Diffusion of labels can be set up as a quadratic function H:

        ![equation](https://latex.codecogs.com/gif.latex?H&space;=&space;\sum_{i}(y_{i}-f{i})^{2}&space;&plus;&space;\\alpha&space;\sum_{i,j}G_{ij}(f_{i}-f_{j})^{2})

        where ```y``` is the initial label vector (```y(i) = 1``` for a protein carrying the
        label, and ```0``` otherwise), and ```f``` is the final label vector post-diffusion.
        The goal is to minimize H by balancing the two counter-acting terms:

        * The first term ```(y-f)``` represents loss of inital label, and we
            want to minimize it (eg the originally-labeled nodes need to retain most of
            their label).

        * The second term ```f(i) - f(j)``` represents the difference in label state
            of connected nodes (```G(i,j)``` is our graph, where ```G(i,j)=1``` if the two nodes
            are connected and ```0``` otherwise). To minimize this term, we want to ensure
            that label information is distributed smoothly over the network. To promote
            smoothing, we can modulate constant ```alpha```.

        Given this linear system, we can solve for ```f``` using the graph diffusion
        kernel. In practise,```f``` is the
        post-diffusion "information content" of the node, and we use it as a measure
        of the node's connectedness to the orignal site of the label.

        ### Kinase network as model system
        ---

        So far, we have explained how to (1) build a network out of
        GO terms, and how to (2) use diffusion to find clusters of connected proteins within the network.

        To test this approach, I chose to model the network of the human kinome. I chose
        to limit things to the kinome for the following reasons:

        * I wasn't sure the browser can render a full network of 20,000 human proteins.
        * Kinases are important to human disease, and are relatively well-studied and annotated

        The network of 327 human kinases (the ones with at least 10 GO terms) is the
        network used in this dashboard. You can diffuse labels between nodes in this
        network to find kinases that are potentially connected.

        ### Does this tool actually predict novel associations? (Validation & Limitations)
        ---
        Some of the kinase network validation can be found in Chapter 3 of [my thesis](https://github.com/ily123/thesis).
        To summarize my findings:

        * Cross-validation across multiple label sets is pretty decent with AUCs > 0.80.
        This shows that the tool is good at modeling existing GO term information as a protein graph.

        * However, prediction of future interactions, as measured by retrospective validation, is not as impressive.
        Retrospective validation produces AUCs for best label cases in the range of 0.6-0.7.
        * Hold-out validation lagging behind cross-validation is a sign of overtraining.

        Why does the GO network overtrain? I think this is because GO term annotations
        are binary. That is, if a protein has no terms linking it to a specific process,
        the model has no data to build connections with. For
        annotation of truly orphan proteins, GO term networks, as presented here,
        are not that useful. For underannotated proteins, it's better to use intrinsic
        properties like primary sequence or structure (something that can also be made
        into networks).

        So in the worst case scenario, you can use this tool as an engaging way to look at
        GO term information. In the best case scenario, you can use it to predict
        novel associations. See the EXAMPLE tab for a fully worked example of how to
        actually use the tool.

        Finally, feel free to take my code and plug in your own networks. If
        I had more time to work on this, that's probably where I'd take this tool.
        """,
    "example_tab": """
        ### Typical usecase
        The goal of the GGID tool is to find clusters of closely connected kinases.

        Let's say you are interested in the protein P53. It's central to several cancers,
        and interacts with many kinases. You read the literature and realize that interaction
        of P53 with its kinases is critical to several cancers.

        In addition to any kinases you already know about,
        you want to find any additional kinases that may interact with P53. So your hypothesis
        becomes:

        ```
        Close network neigbors of P53 kinases are also P53 kinases.
        ```

        Here is how to generate a list of these putative neigbors using my tool:

        ---

        **Step 1: Come up with a list of known P53 kinases.**

        You can do it by reading literature or using a database like [Phosphosite](www.phosphosite.org). In this
        case, lets skip over to [P53's page on Phoshosite](https://www.phosphosite.org/proteinAction.action?id=465&showAllSites=true)
        and then go the *Upstream* tab. This tab lists all upstream regulators of P53.
        One of the categories in the tab is *Kinases, in vitro*. This is the list of
        kinases that phosphorylate P53 in the test tube. Here is the list:
        ```
        ATM, ATR, AurA, AurB, Btk, CDK1, CDK4, CDK5, CDK9, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1, LRRK2,
        MAPKAPK5, NEK2, NuaK1, P38A, P38G, PAK4, PKCD, PLK3, PRPK, SMG1, Src, TAF1, VRK1
        ```

        ---

        **Step 2: Enter data and conduct diffusion experiment:**

        Paste the kinase list into the text box, check cross-validation, and press the
        DIFFUSE button.

        ![alt text](assets/screenshot_2.png)

        ---

        ** Step 3: Analyze initial results:**

        The output of the diffusion experiment is the list of proteins ranked by how
        closely they are connected to the input label set. By default all proteins
        ranked with zscore of 2 or more are plotted on the graph (you can adjust this
        value in the input form).

        So here is what the network of our top hits looks like now:

        ![alt text](assets/screenshot_3a_circle.png)

        The diamonds are the orignal input labels, and circles are the top hits.
        Color of the circle corresponds to the z-score. You can view exact ranks and
        z-scores in the table below the graph (not shown here).

        The circular layout is not super informative, so click on Force-directed layout.
        The result will look something like this:

        ![alt text](assets/screenshot_3b_force_dir.png)

        The graph looks like a bit of a mess, but we begin to see some separation between
        nodes.

        However, for now, let's focus on two things:


        * There are some input labels not connected to anything (boxed in red).
        Why is that? Well, the input set is pretty broad, there are multiple proteins
        so there could be several clusters. The proteins that don't cluster to anything
        are either under annotated, or simply don't connect well to the other labels.


        * The metric of how good the labels connect to each other is the cross-validation
        AUC. You can read about cross-validation in detail [here](https://en.wikipedia.org/wiki/Cross-validation_(statistics)#Leave-one-out_cross-validation).
        We leave one of our input labels unlabeled, and diffuse the remaining N-1 labels
        to see if the left-out label is closely connected to them. We repeat this for
        every label, and summarize results as ROC curve. The close the are under the
        curve, the more connected (more predictive of each other) the input labels are.

        **Why do we care if the labels are connected to each other?** If they are connected
        to each other, then our initial assumption that these labels do the same thing
        (and can therefore predict similar things) is correct. If the labels are not
        connected, the our input set is not informative. It's just a collection of unconnected
        proteins!

        Right now, our AUC is not very high. It's 0.61, meaning the lables are somewhat
        connected, but there is a lot of noise. We can work with the list of the top
        hits, but let's try to fine-tune our inputs.

        ---

        **Step 3b (optional): Adjust the inputs to improve cross-validation performance**

        We have some unconnected proteins in the network. Let's go ahead and get rid
        of them. In this example, those proteins are BTK, CDK9, NEK2, PAK4 and PLK3.
        Remove them from the input set and rerun the experiment.

        New input set:

        ```
        ATM, ATR, AurA, AurB, CDK1, CDK4, CDK5, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1, LRRK2,
        MAPKAPK5, NuaK1, P38A, P38G, PKCD, PRPK, SMG1, Src, TAF1, VRK1

        ```
        Re-run the experiment. The results will look like so:
        ![alt text](assets/screenshot_3c_iteration1.png)

        The basic connectivity of the major cluster hasn't changed, but the isolated
        labels no longer introduce noise into the experiment. (And because we excluded
        those labels from the analysis, the cross-validation performance of the
        remaining labels improved to AUC=0.73.)

        From here, we can continue excising proteins that are not clearly connected
        to a cluster. Two prime targets are SRC and LRRK2. In the image below, I dragged
        them out of the hairball to show how many connections they have:


        ![alt text](assets/screenshot_3c_iteration2.png)

        They have many connections,
        but those connections aren't very specifict. In a synthetic network
        like mine, too many edges suggests that the proteins isn't strongly connected
        to anything. (You can also see it in their low ranks, if you were to look in
        the z-score table, you will find that SRC and LRRK2 are the two lowest
        scoring labels).
        Note: I need to think/read about hub proteins more deeply. So my explanation
        may not be entirely corrent.

        Remove these two from the set, and re-run the experiment. The input label set
        now looks like this:

        ```
        ATM, ATR, AurA, AurB, CDK1, CDK4, CDK5, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1,
        MAPKAPK5, NuaK1, P38A, P38G, PKCD, PRPK, SMG1, TAF1, VRK1
        ```

        The output is now:

        ![alt text](assets/screenshot_3c_iteration3_.png)

        With SRC and LRRK2 gone, we can now see clear clusters in the remaining labels
        and top hits. One last unconnected protein is GRK5. We can remove it, rerun, and
        the final state of the experiment would be this:

        ![alt text](assets/screenshot_3c_iteration4.png)


        The cross-validation AUC is a respectable 0.87, and the hairball is gone. The
        labels and their top hits are now interpretable by a human.

        **Step 5: Interpret the results:**

        Let's circle back to our original hypothesis - clusters of connected proteins
        predict other proteins involved in the same activity by the virtue of being
        connected. We now have a set of input labels that are connected (dimonds, AUC=0.87).

        So then the top hits (proteins that were left unlabled) are our primary targets
        for being involved in some manner with P53. All those of proteins with z-score
        >2 are in the graph.

        You can also see the table at the bottom:

        ![alt text](assets/screenshot_4a.png)

        Note that the labels are present in these tables (initial_state = 1). All those
        proteins with initial_state=0 are the predictions. **To run the experiment without
        LOO (no labels in the output table), uncked LOO box and re-run.**

        You can now pass this list of top hits to your domain/experts experimentalists
        with the following message:

        ```
        Proteins on the top of the list are likely to be associated/interact with P53 or
        P53 related processes: HIPK1, PRKDC, TLK3, BUB1B, etc, in the order of decreasing
        z-score.

        ```
        And you can supply them with the graph, to see if any specific clusters / associations
        jump out.


        **Step 6: Bonus / Loose Ends**

        So what about all those labels we removed? Are they useless? Some of them are
        (in the context of our synthetic GO term network). They are poorly annotated
        and not connected to anything too specific.

        Others are simply involved in some facet of P53 action that's not common to the
        labels that remained in our final set. To get a better idea what those labels
        interact with, you should diffuse them individually.
        """,
}
