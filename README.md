# **Google Drive Link**
https://drive.google.com/drive/folders/1D08TUE-7Jbl2cNmLnid-oX9u14z8g9hm?usp=sharing

# **Off Target Executables**
* Executables for OT are already compiled and located in the Google Drive provided above. Once in the Google Drive, navigate to your OS folder, and executables are located in the "Old" and "New" folders.
* The "Old" folder holds the executables for the old version of OT that uses the original compression format of .cspr files.
* The "New" folder holds the executables for the new version of OT that uses gzip compressed .cspr files, and .db files for repeat data.
* If you want to build the OT source code from scratch you'll need the Boost, Sqlite3, and Zlib libraries. You can follow the following sections for building the libraries and for building OT; otherwise, skip to the "How to use OT" section to see the command line arguments needed for OT.

# **Download and Compile Boost v1.73**

**Mac and Linux**:
1. Download Boost v1.73 .gz folder (in Google Drive)
2. Extract boost files to "boost_1_73_0" folder
3. Prepare the installer by running: `./bootstrap.sh --with-libraries=atomic,chrono,context,coroutine,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,log,math,mpi,program_options,random,regex,serialization,system,test,thread,timer,wave`
4. If on Linux, install prerequisites for the iostreams: `sudo apt-get install libbz2-dev`
5. Build boost with the following command: `sudo ./b2 install`

**Windows (Visual Studio 2017)**:
1. Download Boost v1.73 .zip folder (in Google Drive)
2. Download zlib v1.2.11 .zip folder (in Google Drive)
3. Extract boost files to "boost_1_73_0" folder
4. Extract zlib files to "zlib-1.2.11" folder
5. Open CMD
6. CD to "boost_1_73_0" or wherever the boost files were extracted to
7. Prepare the installer by running: `bootstrap.bat --with-libraries=atomic,chrono,context,coroutine,date_time,exception,filesystem,graph,graph_parallel,iostreams,locale,log,math,mpi,program_options,random,regex,serialization,system,test,thread,timer,wave`
8. Build the boost library with the following command (works with Visual Studio 2017):  `b2 -j15 --toolset=msvc --build-type=complete stage -s ZLIB_SOURCE="C:\Path\To\zlib-1.2.11" -s ZLIB_INCLUDE="C:\Path\To\zlib-1.2.11"`


# **Download And Compile Sqlite3**
**Linux**: 
1. Open terminal
2. Run the following command: `sudo apt-get install sqlite3`

**Mac and Linux (if you did not use apt-get)**:
1. Download sqlite3 source code (in Google Drive)
2. Open terminal
3. CD to sqlite3 source code folder
3. Build using the following command: `gcc -c sqlite3.c`
	* The build command should generate a .o file that will be used with compiling OT.
4. Copy the .o file to the OT source code folder

**Windows (Visual Studio 2017)**:
1. Download sqlite3 source code (in Google Drive)
2. Open "Developer Command Prompt for VS 2017"
3. CD to sqlite3 source code folder
4. Build the library files by running the following command: `lib /DEF:sqlite3.def /OUT:sqlite3.lib /MACHINE:x64`

# **Download and Compile OT**
**Linux (if you used apt-get install sqlite3)**:
1. Download OT source code for Linux (in Repository)
2. Open terminal
3. CD to OT source code folder
3. Run Command to compile OT: `g++ -std=c++11 CasperOffTargetWindows.cpp CSeqTranslate.cpp csprRef.cpp FileOp.cpp gRNA.cpp OffScoring.cpp pch.cpp RefSequences.cpp -pthread -lsqlite3 -lboost_system -lboost_iostreams -o OT`

**Mac and Linux (if you manually built sqlite3 .o file, make sure sqlite3 .o file is in same folder as OT source code)**:
1. Download OT source code for Mac or Linux (in Repository)
2. Open terminal
3. CD to OT source code folder
3. Run command to compile OT: `g++ -std=c++11 CasperOffTargetWindows.cpp CSeqTranslate.cpp csprRef.cpp FileOp.cpp gRNA.cpp OffScoring.cpp pch.cpp RefSequences.cpp -pthread sqlite3.o -lboost_system -lboost_iostreams -o OT`

**Windows (Visual Studio 2017)**:
1. Download OT source for Windows (in Repository)
2. Open Visual Studio 2017
3. Create a New Project
4. Import the OT source code files
5. Go to Project->Properties and do the following:
	* Set Configuration to "All Configurations"
	* Set Platform to "All Platforms"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\boost_1_73_0;C:\Path\To\sqlite3;` to "Include Directories"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\boost_1_73_0\stage\lib;C:\Path\To\sqlite;` to "Include Libraries"
	* In C/C++->General add `C:\Path\To\boost_1_73_0;C:\Path\To\sqlite3;` to "Additional Include Directories" and "Additional #using Directories"
	* In Linker->General add `C:\Path\To\boost_1_73_0\stage\lib;C:\Path\To\sqlite;` to "Additional Library Directories"
	* In Linker->Input add `C:\Path\To\sqlite\sqlite3.lib;` in "Additonaly Dependencies
	* Make sure when you add these paths that there are ';' seperating all paths/object in the line.
6. For debugging, make sure Debug and x86 are selected before running
	* If debugging, make sure you set the debugging command arguments. See "How to run OT" below.
7. For compiling, make sure Release and x86 are selected before running
	

# **How to run OT**
* CD to the directory containing the OT executable
	* The command line arguments for the new OT are as follows: `query_file cspr_file db_file output_file CASPERinfo_file num_mismatches threshold detailed_output_bool avg_output_bool`
	* The command line arguments for the old OT are as follows: `query_file is_compressed_bool cspr_file output_file CASPERinfo_file num_mismatches threshold detailed_output_bool avg_output_bool`
* Example command (New OT): `./OT query.txt ath_asCas12.cspr ath_asCas12_repeats.db output.txt CASPERinfo 5 0.05 True False`
* Example command (Old OT): `./OT query.txt True ath_asCas12.cspr output.txt CASPERinfo 5 0.05 True False`
