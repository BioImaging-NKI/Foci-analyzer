#@ Integer smoothTracesStep (label = "Smooth traces with (radius)", value=2, min=0);

Table.create("Foci counts over time");
Plot.create("Foci counts", "time (frames)", "foci count");

cellIDData = Table.getColumn("Cell ID", "Results");
focicountData = Table.getColumn("Foci count ch3", "Results");
timeData =  Table.getColumn("filename", "Results");

cellIDs = uniqueValuesInArray(cellIDData);
uniqueNrCells = cellIDs.length;

setBatchMode(true);

//Get the Glasbey LUT from a random image
newImage("Untitled", "8-bit black", 32, 32, 1);
run("glasbey_on_dark");
getLut(reds, greens, blues);
close();

for (i = 0; i < uniqueNrCells; i++) {
	fociTrace = lookup(cellIDData, i, focicountData);
	timeTraceString = lookup(cellIDData, i, timeData);
	timeTrace = newArray(timeTraceString.length);
	for (k = 0; k < timeTraceString.length; k++) {
		timeTrace[k] = substring(timeTraceString[k], indexOf(timeTraceString[k], "__t=")+4,indexOf(timeTraceString[k], "__t=")+7);
	}
	Table.setColumn("Cell "+cellIDs[i], fociTrace);
	color = getLabelColor(i, cellIDs.length);
	Plot.setColor(color);
	fociTraceSmoothed = smoothArray(fociTrace, smoothTracesStep);
	Plot.add("line", timeTrace, fociTraceSmoothed);
}
Table.update;
Plot.show();
Plot.setLimits(0, NaN, NaN, NaN)
setBatchMode("show");




//Returns, as array, all values in the array return_array at the same indices as lookup_value in the ref_array
function lookup(ref_array, lookup_value, return_array) {
	indices = indexOfArray(ref_array, lookup_value);
	return_values = newArray(indices.length);
	for (i=0; i<indices.length; i++) {
		return_values[i] = return_array[indices[i]];
	}
	return return_values;
}


//Returns, as array, the indices at which a value occurs within an array
function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
	return newArray(0);
}


//Returns an array with unique values in the array
function uniqueValuesInArray(array) {
	count=0;
	uniqueArray = newArray(array.length);
	for (a=0; a<lengthOf(array); a++) {
		if(!occursInArray(Array.trim(uniqueArray, count), array[a]) && !matches(array[a],"None")) {
			uniqueArray[count]=array[a];
			count++;
		}
	}
	return Array.trim(uniqueArray, count);
}


//Returns true if the value occurs within the array
function occursInArray(array, value) {
	for(i=0; i<array.length; i++) {
		if(array[i] == value) return true;
	}
	return false;
}


function getLabelColor(label, nrOfLabels) {
	if(nrOfLabels >= 255) {
		color1 = IJ.pad(toHex(reds[label/nrOfLabels*255]),2);
		color2 = IJ.pad(toHex(greens[label/nrOfLabels*255]),2);
		color3 = IJ.pad(toHex(blues[label/nrOfLabels*255]),2);
	}
	else {
		color1 = IJ.pad(toHex(reds[label]),2);
		color2 = IJ.pad(toHex(greens[label]),2);
		color3 = IJ.pad(toHex(blues[label]),2);		
	}
	labelColor = "#"+color1+color2+color3;
	return labelColor;
}


//Returns the moving average smoothed array with a radius of 'step'. Edges are handled symmetrically (less data, so more noise towards the edges)
function smoothArray(array, step) {
	if(step>0) {
		smoothed_array = newArray(array.length + step + 1);
		for(i=0; i<=array.length + step; i++) {
			if(i < 2*step) 				piece = Array.slice(array, 0, i+1);
			else if (i < array.length)	piece = Array.slice(array,i-step-1,i+step-1);
			else 						piece = Array.slice(array,i-step-1,array.length+1);
			Array.getStatistics(piece, min, max, mean, stdDev);
			smoothed_array[i] = mean;
		}
		smoothed_array = Array.slice(smoothed_array, round(step/2), array.length + round(step/2));
		return smoothed_array;
	}
	else return array;
}
