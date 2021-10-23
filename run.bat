:: Environment variables

set input-file=tests\PFEM3D_Metafor\damBreak\Wcompressible\input_fsi.py
call "C:\Local\LibsVs2017Py3\LibsVS2017Py3.cmd"

:: Clean output folder

rd /s /q workspace
mkdir workspace

:: Runs the code

python run.py %input-file%