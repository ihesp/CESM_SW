       We have done some efforts on refactoring and optimizing the Community
       Earth System Model(CESM) on the Sunway TaihuLight supercomputer, which
       uses a many-core processor that consists of management processing
       elements (MPEs) and clusters of computing processing elements (CPEs).
       To map the large code base of CESM.
       For general questions about those codes, please contact
       [Pilot National Laboratory for Marine Science and Technology (Qingdao)]
       (http://www.qnlm.ac/) or [National supercomputing Center in Wuxi](http://www.nsccwx.cn/)

       The Community Earth System Model version 1.3 (CESM1.3)

                     http://www2.cesm.ucar.edu

       See the CESM web site for documentation and information.


For commits to the CESM svn repo

- 1) check out the latest ccsm4 tag but do not populate external directories
```shell
> svn co --ignore-externals $SVNREPO/cesm1/tags/cesm1_3_alpha## my_sandbox
> svn co -N $SVNREPO/cesm1/tags/cesm1_3_alpha## my_sandbox
```

- 2) modify the file SVN_EXTERNAL_DIRECTORIES with any changes to
component tags
```
> cd  my_sandbox
> emacs SVN_EXTERNAL_DIRECTORIES
```

- 3) set the property for the external definitions - don't forget the dot at the end
```shell
> svn propset  svn:externals  -F SVN_EXTERNAL_DIRECTORIES .
```

- 4) populate your sandbox with the external files
```shell
> svn update
```

- 5) test

- 6) document your mods
```shell
> cp  ChangeLog_template  ChangeLog.new
> cat  ChangeLog  >>  ChangeLog.new
> emacs  ChangeLog.new
> mv ChangeLog.new  ChangeLog
```

- 7) copy your sandbox to a new tag in the repository
```shell
> svn copy . $SVNREPO/ccsm4/tags/cesm1_3_alpha## -m "created tag cesm1_3_alpha##"
```
