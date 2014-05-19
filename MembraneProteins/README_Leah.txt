5/19/14

Hey Leah,
This is the READ ME to generate your half of the membrane protein sector database.

First and foremost, sync the github repository to make sure you have the latest version of all our files.

Open matlab and add to the path the whole github repository and all subfolders:
	- in the Matlab file tree (left panel), go one level up from the github repository
	- right click on the github folder
	- click add to path > folder and subfolders
Now you should be able to access any function stored in the github repository.

Still on the matlab file tree, change your working directory to the following directory: github/MembraneProteins.
You should now see all the files in MembraneProteins in the file tree on the left panel.

Double click on membraneSectorDBScript_Leah.m to open the script that will create the sector database. I have modified your script to run on your half of the database and create variables with you name appended.

Run the script. Hopefully, it will work. If it doesn't and you can't troubleshoot it yourself, let me know ! 
The script will keep all the fasta files and store them in the MembraneProteins folder, don't erase them ! We might need them later.

If you could try to run this overnight, it will be great to have it by tomorrow.

