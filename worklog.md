#### jun 28 session
    Wrote validator for id input in the app.
    App now filters out anything that is not a valid kinase HUGO id.

#### jun 29 session
    Modified app to feed diffusion script actual string IDs.
    Hooked up diffusion script with id coverter.

#### jun 30 session
    App now shows actual diffusion output as table.
    Scores now include z-score, and hide the initial y=1 proteins (labels).
    Move up submit bar to top of screen, so table shows up in middle of screen when you diffuse.
    Added title to tab, added favicon.

### jul 1 session
    Created SimilarityMatrix container class that holds both the network and the protein ids.
    Inputs now correctly map to the network.
    Formatted Data table (pagination and precision) -- almost good, Dash seems to hide trailing 0s


### jul 4 session
    Wrote basic LOO code. It works but is ugly.

### jul 8 session
    trying to get p53 loo to work // fixing bugs
    loo auc for p53 kinases is strangely low (0.60), something is broken

### jul 11 session
    filtered out proteins with low annotation number (auc still ~0.65)
    examined BMA network derived by ffsim - edge disto a lot smoother
    it appears that using max(sim(a,b)) edges creates lots of equal
        edges that don't threshold (numpy partition does not break ties),
        and it desentitizes the network
    implemented BMA - YAY! Network  works now, AUC = 0.87!

### jul 15 session
    made png figure of post-diffusion distro of z-scores
    inserted figure into web-app /ugly/ - i don't understand html/div positioning logic
    replaced png figure with dash's own 'graph' component

### jul 19 session
    fixed DF output to presented data sorted
    added zscore and rank to web app output
    added .css style sheet
    added test outline to each div
    inserted cytoscape network view into web-app

### jul 20 session
    rewrite tracking protocol
    construct network view for label and best results
