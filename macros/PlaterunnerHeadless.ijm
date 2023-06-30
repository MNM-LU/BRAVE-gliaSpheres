macro "AdjustPlaterunnerNeuronTau" {
run("Image Sequence...", "open=[/home/rstudio/seqFiles/Images/20170418_Tau_Round3LE/ex02em02/001_R3LE-Plate-01_A01_GFP.tif] sort");
    args = split(getArgument(),":"); 
    print(args[0]);
    print(args[1]); 
    run("Image Sequence...", "open=[/home/rstudio/seqFiles/Images/20170418_Tau_Round3LE/ex04em05/001_R3LE-Plate-01_A01_mCherry.tif] sort");
    open("/home/rstudio/seqFiles/Images/20170418_Tau_Round3LE/Transformations.xls");
	selectWindow("ex02em02");
    run("Subtract Background...", "rolling=300 stack");
	headings = split(String.getResultsHeadings); 
	for (row=0; row<nResults; row++) { 
		line = "";
		setSlice(getResult(headings[0],row));
		moveX= 142-getResult(headings[1],row);
		moveY= 142-getResult(headings[2],row);
		run("Translate...", "x="+ moveX + " y="+ moveY + " interpolation=None slice");
	} 

	selectWindow("ex04em05");
    run("Subtract Background...", "rolling=300 stack");
	headings = split(String.getResultsHeadings); 
	for (row=0; row<nResults; row++) { 
		line = "";
		setSlice(getResult(headings[0],row));
		moveX= 142-getResult(headings[1],row);
		moveY= 142-getResult(headings[2],row);
		run("Translate...", "x="+ moveX + " y="+ moveY + " interpolation=None slice");
	} 

	run("Clear Results");
	selectWindow("ex02em02");
	run("Set Measurements...", "area mean min redirect=None decimal=3");
	saveSettings;
	makeOval(142, 142, 1764, 1764);
	for (n=1; n<=nSlices; n++) {
          setSlice(n);
          run("Measure");
        }
    	restoreSettings;
    	saveAs("Results", "/home/rstudio/output/NeuronTau_inVitro_GFP_Results.xls");
    	run("Clear Results");

	selectWindow("ex04em05");
	run("Set Measurements...", "area mean min redirect=None decimal=3");
	saveSettings;
	makeOval(142, 142, 1764, 1764);
	for (n=1; n<=nSlices; n++) {
          setSlice(n);
          run("Measure");
        }
    	restoreSettings;
    	saveAs("Results", "/home/rstudio/output/NeuronTau_inVitro_mCherry_Results.xls");
    	run("Clear Results");


	selectWindow("ex02em02");
	makeOval(142, 142, 1764, 1764);
 	run("Window/Level...");
	setMinAndMax(0, 4000);
	run("Apply LUT", "stack");
	run("Close");
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");

	selectWindow("ex04em05");
	makeOval(142, 142, 1764, 1764);
 	run("Window/Level...");
	setMinAndMax(0, 6000);
	run("Apply LUT", "stack");
	run("Close");
	run("Make Inverse");
	setForegroundColor(0, 0, 0);
	run("Fill", "stack");

	run("Merge Channels...", "c1=ex04em05 c2=ex02em02");
	makeOval(142, 142, 1764, 1764);
	saveAs("Tiff", "/home/rstudio/output/Round2_GFP-mCherry.tif");
    run("Quit");
}

