/*
 * Macro to merge results files.
 * Bram van den Broek, Netherlands Cancer Institute, 2013-2021
 * 
 */

#@ String(label = "String that ends the filenames, including '.' (will be cut off in the resulting name)", value="__Foci_results.tsv") filestring
#@ File (label = "Results folder", style = "directory") dir
#@ String (label = "Output file name", value="Results_all_files") outputfilename


dir = dir+File.separator;

print("\n--- Combining '"+filestring+"' files ---");
run("Clear Results");
run("Input/Output...", "save_column");

//Find all result files in directory
//dir = getDirectory("Choose a directory containing the experiment.");
showStatus("reading directory");
file_list = getFileList(dir); //get filenames of directory
showStatus("");

print(dir);

//make a list of tables (files that end on 'string'.
j=0;
table_list=newArray(file_list.length);
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],filestring) && !matches(file_list[i], outputfilename) && !matches(file_list[i], "All_Foci_Results_per_cell")  && !matches(file_list[i], "All_Foci_Statistics")) {
		table_list[j] = file_list[i];
		j++;
	}
}
table_list = Array.trim(table_list, j);	//Trimming the array of images
if(table_list.length==0) exit("No results files found");
else print(table_list.length+" results files found in "+dir+":");

//Get column headers from the first file in the list
text_file = File.openAsString(dir+table_list[0]);
lines = split(text_file,"\n");
if(substring(filestring, lastIndexOf(filestring, ".")) == ".csv") headers = split(lines[0],",");
else headers = split(lines[0],"\t");

//import all results into the results window and save
j=1;
m=0;
for(i=0; i<table_list.length; i++) {
	print(table_list[i]);
	text_file = File.openAsString(dir+table_list[i]);
	lines = split(text_file,"\n");
	for(j=1;j<lines.length;j++) {
		if(substring(filestring, lastIndexOf(filestring, ".")) == ".csv") values = split(lines[j],",");
		else values = split(lines[j],"\t");
//		if(j==1) setResult("filename",m,substring(table_list[i], 0, lengthOf(table_list[i])-12));	//Only print file name at the first occurrance
//		else setResult("filename",m,"");
		setResult("filename",m,substring(table_list[i], 0, lengthOf(table_list[i])-lengthOf(filestring)));
		for(k=0;k<headers.length;k++) {
			if(parseFloat(values[k]) == values[k]) setResult(headers[k],m,parseFloat(values[k]));
			else setResult(headers[k],m,values[k]);
		}
		m++;
	}
	updateResults();
}
selectWindow("Results");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
saveAs("text", dir+outputfilename+"__(t="+IJ.pad(hour,2)+"_"+IJ.pad(minute,2)+"_"+IJ.pad(second,2)+").tsv");
print("\nResults File saved as "+dir+outputfilename+"__(t="+IJ.pad(hour,2)+"_"+IJ.pad(minute,2)+"_"+IJ.pad(second,2)+").tsv");
