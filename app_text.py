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
        In case of disease, this allows you to expand the candidate list of theraupeutic targets, and arrive at
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
        with labels such as "protein kinase", "signaling", "cell wall",
        and so on. These labels are  known as "GO terms" and they describe a protein's
        job within the cell. Because proteins that do
        similar jobs have similar GO terms, GO term similarity can be readily
        converted into a network/graph (which we can then learn from).

        To build a protein network from GO terms, we follow these steps:

        1. Given a set of proteins, we measure GO term similarity between each possible pair
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
            general high-level term, and one of its children is the more specific ```negative regulation of signaling```.
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

            In practice, the most informative common ancestor is usually the lowest common ancestor in the tree.
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
            partner proteins as edges, and discard all other pairs. This produced a protein
            graph where each protein is connected to 5 proteins that it is most similar to.

        3. To learn from this network, we then use information diffusion
        as described below.

        ### Graph Information Diffusion
        ---

        To extract information from the network, we can simply examine it by hand. Given a protein, what
        are its closest neighbors? This works, but in practice gets out of hand quickly
        when you want to examine more than a single protein at the same time, or when the number of edges
        is too high. So it's best to use an algorithmic approach, such as Information
        Diffusion.

        Intuitively, Information Diffusion is exactly what it sounds like. We label
        some nodes in the network with information, then let information diffuse to other
        nodes using edges as conduits. At the end of the diffusion experiment, nodes that accrue the
        most information are the ones that are most
        closely connected to the originally-labeled nodes.

        Using this method, we can discover clusters of proteins that do the same job.
        For example, label all proteins associated with ```disease X```, diffuse
        that label, and measure amount of information in each node post diffusion.
        Rank the originally unlabled nodes by the amount of information they retained.
        Proteins at the top of the ranked list are likely to be involved in
        ```disease X```.

        For formal mathematic definition of diffusion see the methods section
        of [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3001439/).

        ---
        However, here is a quick summary from someone who isn't that great at linear
        algebra:

        Diffusion of labels can be set up as a quadratic function H:

        ![equation](https://latex.codecogs.com/gif.latex?H&space;=&space;\sum_{i}(y_{i}-f{i})^{2}&space;&plus;&space;\\alpha&space;\sum_{i,j}G_{ij}(f_{i}-f_{j})^{2})

        where ```y``` is the initial label vector (```y(i) = 1``` for a protein carrying the
        label, and ```0``` otherwise), and ```f``` is the final label vector post-diffusion.
        The goal is to minimize H by balancing the two counter-acting terms:

        * The first term ```(y-f)``` represents loss of initial label, and we
            want to minimize it (eg the originally-labeled nodes need to retain most of
            their information).

        * The second term ```f(i) - f(j)``` represents the difference in label state
            of connected nodes (```G(i,j)``` is our graph, where ```G(i,j)=1``` if the two nodes
            are connected and ```0``` otherwise). To minimize this term, we want to ensure
            that label information is distributed smoothly over the network. To promote
            smoothing, we can modulate constant ```alpha```.

        Given this linear system, we can solve for ```f``` using the graph diffusion
        kernel. In practice, ```f``` is the
        post-diffusion "information content/label/state" of the node, and we use it as a measure
        of the node's connectedness to the original site of the label. In the final tally of
        ```f``` scores, we also convert them to z-scores for convenience/generalizability.

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
        The goal of the GGID tool is to find clusters of closely connected kinases
        within the GO-term similarity network of the human kinome.

        When would you want to do that? Here is a typical scenario:

        Let's say you are interested in the protein P53, because it's involved in a
        cancer you care about. You read the literature on P53 and realize that it
        interacts with many kinases, and these interactions are critical to disease
        progression. In order to cure the P53-related cancer, you decide to test whether any of
        these kinase interactions can be used as drug targets.

        However, before you embark on your drug screen, you wonder if your list of P53
        kinases (the one you got from the literature) is complete. Are there any
        other putative P53 kinases not known to science that you could test?

        To answer this question, you can use the GGID tool!
        The basic hypothesis the tool provides for you is:

        ```
        Within the network of human kinases, neighbors of P53 kinases are also P53 kinases.
        ```

        Given a set of input kinases (in this case, known P53 kinases), the tool
        labels these kinases in the network, and then diffuses the label to other nodes in the
        network. The tool then measures how much information diffused to each
        previously-unlabeled node, and uses it as proxy to judge how well each node is
        connected to the initial label set (see THEORY tab).
        Kinases at the top of the post-diffusion label ranked list are,
        hypothetically, connected to P53. So after verifying,
        you can include them in your drug screen.

        ---

        Before I continue with the P53 example, I'd like to emphasize
        that the tool is completely generic. In this example, we are using P53 kinases,
        but input kinase labels can be any set of kinases that you think share
        some relavant property (kinases of some other protein, kinases involved in a specific
        signaling pathaway, kinases involved in a specific disease, etc).

        ### How to use the tool (to suggest novel P53 kinases)

        **Step 1: Come up with a list of known P53 kinases.**

        Before you start, you need to identify the list of known positives, i.e. known
        P53 kinases.

        You can do it by reading the literature or using a database like [Phosphosite](www.phosphosite.org).
        For the purposes of this demo, let's skip over to [P53's page on Phosphosite](https://www.phosphosite.org/proteinAction.action?id=465&showAllSites=true)
        and then go the *Upstream* tab. This tab lists all upstream regulators of P53.
        One of the categories in the tab is *Kinases, in vitro*. This is the list of
        kinases that are known to phosphorylate P53 in the test tube. Here is the list:
        ```
        ATM, ATR, AurA, AurB, Btk, CDK1, CDK4, CDK5, CDK9, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1, LRRK2,
        MAPKAPK5, NEK2, NuaK1, P38A, P38G, PAK4, PKCD, PLK3, PRPK, SMG1, Src, TAF1, VRK1
        ```

        ---

        **Step 2: Enter known kinases and conduct diffusion experiment (with cross-validation):**

        Paste the kinase list into the text box, check cross-validation, and press the
        DIFFUSE button.

        ![example pic](assets/example/example_pic1.png)

        The cross-validation option is important, and should always be checked if your
        input set is larger than 1 protein. It conducts a cross-validation experiment
        on your input label set. We will discuss cross-validation in more detail below.

        ---

        ** Step 3: Analyze initial results:**

        The output of the diffusion experiment is the list of proteins ranked by how
        closely they are connected to the input label set (measured as label z-score, see
        THEORY for details). By default all proteins with z-score >= 2 are plotted on the graph.

        So here is what the graph of our top hits looks like now:

        ![example pic](assets/example/resized-example_pic2.png)

        Diamonds are the orignal input labels, and circles are the unlabeled proteins with
        post-diffusion z-score of 2 or more. Color
        of the node corresponds to their rank/z-score. You can view z-score of each
        protein by clicking on it. For example, CDK2:

        ![example pic](assets/example/resized-example_pic3.png)

        Meanwhile, the circular layout is not super informative, so click on Force-directed layout.
        The result will look something like this:

        ![example pic](assets/example/resized-example_pic4.png)

        The graph looks like a bit of a mess, but we begin to see some separation between
        nodes.

        For now, let's focus on two things:


        * There are some input labels, boxed in red, not connected to the main cluster.
        Because they are not connected to the rest of the input set, information
        does not diffuse to them so they score low in cross-validation.
        Especially low are SRC, PAK4 (tied for rank=313), and BTK (rank=306).

         Why are these proteins disconnected? It's possible that they are underannoated,
         and and thus don't have strong connections to the more well-annotated proteins in the set.
         Alternatively, they can be part of their own functional cluster. The input set
         of P53 proteins is quite large, and it's likely that there are multiple
         subclusters that do different jobs.

        * The metric for how well the input labels connect to each other is the cross-validation
        AUC. You can read about cross-validation in detail [here](https://en.wikipedia.org/wiki/Cross-validation_(statistics)#Leave-one-out_cross-validation).

         To conduct cross-validation, we leave one of the input labels unlabeled, and diffuse the remaining N-1 labels
        to see if the left-out label is closely connected to the set. If it is, it will
        get a high post-diffusion score. We repeat this for
        every label, and summarize results as a ROC curve. The closer the area under
        the curve is to 1, the better. You can see the curve by clicking on "Cross-validation AUC"
        collapsible.

        **Why do we care if the labels are connected to each other?** If they are connected
        to each other, then our initial assumption that these labels do the same thing
        (and can therefore predict other proteins that are involved in the process) is correct. If the labels are not
        connected, the our input set is not informative. It's just a collection of unconnected
        proteins! Cross-validation allows us to check.

        Right now, AUC is 0.77, meaning the labels are somewhat
        connected, but there is some noise. Let's try to fine-tune our inputs.

        ---

        **Step 4: Adjust the inputs to improve cross-validation performance**

        We have 3 unconnected proteins in the network (SRC, PAK4, and BTK). Let's remove
        them from the input set and rerun the experiment. The new input set is below:

        ```
        ATM, ATR, AurA, AurB, CDK1, CDK4, CDK5, CDK9, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1, LRRK2,
        MAPKAPK5, NEK2, NuaK1, P38A, P38G, PKCD, PLK3, PRPK, SMG1, TAF1, VRK1
        ```
        The results will look like so:
        ![example pic](assets/example/resized-example_pic5.png)

        Because the isolated labels no longer introduce noise into the experiment,
        the AUC is now at a respectable 0.91. Additionally, there is now a visible
        connection between the GRK5-DAPK1 pair to the rest of the cluster via MAP3K12, which
        was previously hidden because MAP3K12 had a z-score slightly below 2.


        From here, we can continue excising proteins that are not strongly connected
        to the cluster. Two prime targets are NEK2 and CDK9 (both rank outside top 100).
        The updated input set is below

        ```
        ATM, ATR, AurA, AurB, CDK1, CDK4, CDK5, Chk1, Chk2, CK1A, CK2B,
        DAPK1, DNAPK, DYRK1A, DYRK2, ERK1, ERK2, GRK5, HIPK2, JNK1, JNK2, LKB1,
        LRRK2, MAPKAPK5, NuaK1, P38A, P38G, PKCD, PLK3, PRPK, SMG1, TAF1, VRK1
        ```
        and the result looks like so:

        ![example pic](assets/example/resized-example_pic7.png)

        The AUC is now 0.95, and we see good sub-clustering of input labels and top hits.
        At this point, we should stop optimizing, and interpret the results.

        **Step 5: Interpret the results**

        Let's circle back to our original hypothesis - network connection predicts
        functional similarity. We now have a solid set of input labels that are strongly
        connected (diamonds, AUC=0.95). By diffusing them, we get a list of other kinases closely connected to
        them in the network. These are the kinases in the graph view.


        In addition to the graph view, you can see the z-score for every kinase
        in the table under "Full z-score and ranks table":

        ![example pic](assets/example/resized-example_pic8.png)

        Note that the input labels (boxed in red) are present in the results table (initial_state = 1), if
        you are running cross-validation. Once you are happy with your input set,
        you can uncheck cross-validation, and the table will only feature unlabeled
        proteins.

        The top-scoring kinases in the table likely functionally connected to P53.
        This is the list of putative P53 kinases. Ideally, you can pass it to your wet-lab
        collaborators and ask them to test whether any of these proteins actually
        interact with P53 in-vitro.
        """,
}
