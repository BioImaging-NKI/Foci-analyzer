/*
 * Macro to merge results files.
 * Bram van den Broek, Netherlands Cancer Institute, 2013-2021
 * 
 */

#@ String(label = "String that ends the filenames, including '.' (will be cut off in the resulting name)", value="_ch1.tsv") fileString
#@ File (label = "Results folder", style = "directory") dir
dir = dir+File.separator;

print("\\Clear");
run("Clear Results");
run("Input/Output...", "save_column");

//Find all result files in directory
//dir = getDirectory("Choose a directory containing the experiment.");
showStatus("reading directory");
file_list = getFileList(dir); //get filenames of directory
showStatus("");

//make a list of tables (files that end on 'string'.
j=0;
table_list=newArray(file_list.length);
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],fileString) && !matches(file_list[i], "Results_all_files.tsv")) {
		table_list[j] = file_list[i];
		j++;
	}
}
table_list = Array.trim(table_list, j);	//Trimming the array of images
if(table_list.length==0) exit("No results files found");
else print("Experiment contains "+table_list.length+" result files:");

//Get column headers from the first file in the list
text_file = File.openAsString(dir+table_list[0]);
lines = split(text_file,"\n");
if(substring(fileString, lastIndexOf(fileString, ".")) == ".csv") headers = split(lines[0],",");
else headers = split(lines[0],"\t");

//import all results into the results window and save
j=1;
m=0;
for(i=0; i<table_list.length; i++) {
	print("appending file: "+table_list[i]);
	text_file = File.openAsString(dir+table_list[i]);
	lines = split(text_file,"\n");
	for(j=1;j<lines.length;j++) {
		if(substring(fileString, lastIndexOf(fileString, ".")) == ".csv") values = split(lines[j],",");
		else values = split(lines[j],"\t");
//		if(j==1) setResult("filename",m,substring(table_list[i], 0, lengthOf(table_list[i])-12));	//Only print file name at the first occurrance
//		else setResult("filename",m,"");
		setResult("filename",m,substring(table_list[i], 0, lengthOf(table_list[i])-lengthOf(fileString)));
		for(k=0;k<headers.length;k++) {
			if(parseFloat(values[k]) == values[k]) setResult(headers[k],m,parseFloat(values[k]));
			else setResult(headers[k],m,values[k]);
		}
		m++;
	}
	updateResults();
}
selectWindow("Results");
saveAs("text", dir+"Results_all_files.tsv");
print("\nAppended Results File saved as "+dir+"results_appended.tsv");
