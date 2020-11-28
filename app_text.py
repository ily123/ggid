"""Module to keep all of the app's text-heavy elements."""
# why not keep the text in a csv file or database?


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
    "example": """
        Let's say you'd like to study protein P53. It's central to several cancers,
        and interacts with many kinases. You have a list of human kinases that
        interact with P53, and you would like to know if there are any more.

        Step 0: Find existing kinases by looking in the Phosposite database
        Step 1: Enter these kinases in teh form
        Step 3: Press diffuse

        Recommended:

        Step 4

        """,
}
