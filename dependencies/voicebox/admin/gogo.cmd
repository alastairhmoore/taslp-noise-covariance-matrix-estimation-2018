set htmpath=C:\work\www\dmb\voicebox

echo creating zip file

del %htmpath%\voicebox.zip
"C:\Program Files\7-Zip\7z.exe" a -tzip %htmpath%\voicebox.zip ..\*.m ..\*.mat ..\*.exe

echo copying matlab documentation

copy ..\doc\*.* %htmpath%\doc
copy ..\doc\voicebox\*.* %htmpath%\doc\voicebox

del ..\*.asv
del ..\*.bak

