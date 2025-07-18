@echo off
set "ff=%~1"
set "prot=%~2"
set "steps=%3"
set "obabel=.\heya\Library\bin\obabel.exe"
set "obminimize=.\heya\Library\bin\obminimize.exe"
set "BABEL_DATADIR=.\heya\share\openbabel"

echo "using %ff% forcefield"
for %%f in (*.sdf) do (
    echo Running command 1 for %%f...
    obabel "%%f" -O "%%~nf_merged.sdf" -p %prot%  
    
    echo Running command 2 for %%f...
    obminimize -o mol2 -sd -ff %ff% -n %steps% -h "%%~nf_merged.sdf"  > "%%~nf_minimized.mol2" 
    
    echo Running command 3 for %%f...
    obabel "%%~nf_minimized.mol2" -O "%%~nf.pdbqt" -m 

    
    move ".\*.pdbqt" ".\output" >nul

)

del  ".\*.sdf" 2>nul
del ".\*.mol2" 2>nul