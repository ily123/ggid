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
        ### The problem
        ---
        Proteins don't work in isolation. They form networks and work in tandem.

        This means that once you've identified
        a protein associated with some process (ex: a disease), you can query that
        protein's network to find additional proteins that are involved.
        Once you understand the network and the role of all relevant protein players,
        you can design effective theurapeutics.

        Construction of such networks, and learning from them, is one of the
        subfields in computational biology. The industry standard for collecting and
        visualizing these data is the [STRING](https://string-db.org/) database. There
        you can find all sorts of protein networks build on top of various data (physical
        interactions, co-expression of proteins, co-expression of protein mRNA, etc).

        ### GO term similarity network
        ---

        In this dashboard, I created a novel way of presenting another type of
        information: GO terms. GO terms come from the [Gene Ontology](http://geneontology.org/)
        database. They are tags assigned to proteins in order to describe proteins job within the cell (very similar to how
        IMDB assigns tags to movies). Proteins that do similar jobs have similar GO terms. We can build a network from that!

        To build a network, we compare proteins to each other, in an all-vs-all manner, to
        produce a similarity matrix that show how similar any pair of proteins are. We then
        reduce the similarity matrix to a graph (network) by dropping out connections
        that are weak (ie the two proteins are not similar). What we end up with is a network.

        The circle pic you see loaded in the first tab is the GO term network for 350 human
        kinases, build exactly like I described. (For those who want the nitty gritty details:
        similarity between two proteins was measured using Resnik similarity, the network
        was constructed by dropping out all edges outside of the top 10. The proteins, beforehand,
        were sanitized by removing those with fewer than 10 annotations.)


        ### Diffusion
        ---

        To extract information from this network, we can simply examine it by hand. Given a protein, what
        are its closest neighbors? That works, but gets out of hand really quickly when you
        want to look at more than just 1 protein. (Let's say you have 10 proteins involved
        in some disease, and you want to see all of their neighbors.)

        To help us, there are a bunch of algorithms that
        simplify the process of finding neighors within networks. The one I am using
        is called Information Diffusion. It was originally developed in the 1990s to
        model diffusion of heat through metal, but has since been deployed to protein
        networks as well.

        The "heat" in our case is information (a binary label: 1 if a protein does X, and 0 otherwise).
        We diffuse the label over the network of proteins. Post-diffusion, the proteins
        that have the highest label content are the ones most closely connected to the
        source.

        ### Kinases
        ---

        Ok, that's all great. But this app only includes a network of kinases.
        why did you build a kinase network specifically?

        I chose kinases for this demo app for technical reasons:

        * I wasn't sure the browser can render a full network of 20,000 human proteins
        * Meanwhile, the number of kinases is realtively small (just 500),
        so they are easy to model
        * Plus, kinases are very important in human disease (they drive signaling
        by phosphorylating other proteins), and as a result well-studied and annotated

        You can apply the framework described here to the entire human proteome (or other actors
        like RNAs).


        ### Does this tool actually work for predicting novel associations?
        ---

        This tool is good at presenting networks of existing information (ie cross-validation
        is pretty decent >= 0.80). That being said, prediction of future interactions
        from present data is pretty iffy (retrospective auc ~0.6-0.7). I am not showing any data here,
        just recalling the experimental results from the waning days of my training.
        I think some of it is in my thesis.

        As to why this model doesn't generalize well (i.e. it overtrains)  - I think this is
        because of the binary nature of GO term annotations. If two proteins are found
        to be associated, but their GO terms are completely different before the discovery
        is made, the model isn't going to be aware the relationship. The model predicts two proteins
        that already share GO terms as functionally related. If there is no similarit to
        cease on, the model simily won't create an edge in th network.

        ALL THIS NEEDS EDITING // FIRST DRAFT


        **If the tool can't predict the future, wth did you make this dashboard??**

        I wanted to write a network app. I have multiple ideas for using the
        basic approach in other domains. And now I have an app for that. Just
        have to change the underlying network and edit some text, and I am good to go!
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
