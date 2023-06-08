# README with useful Git commands

## Cloning lineqGPR and creating a new branch:
1. git clone https://github.com/anfelopera/lineqGPR
2. To create a new branch with name "dev-lastname": git checkout -b dev-lastname
    - To see the existing branches: git branch -a
3. To change the working branch: git checkout dev-lastname

**Note.** to remove a branch: git branch --delete <branchname>


## Pulling updates from the remote/origin/master to master:
1. git fetch origin
2. git checkout master
3. (a) git pull

**Note:** to force the pull from master, instead of (3a), you may use: 3b. git pull https://github.com/anfelopera/lineqGPR master --force

## Merging your branch from master:
1. git checkout dev-lastname
2. (a) git merge master
3a. git push (for the first time: git push --set-upstream origin dev-lastname)

**Note:** to recover the branch from master without tweaking it, instead of (2a), you may use: 2b. git rebase master

## Checking the status of the branch:
- git status

## Removing unsaved changes:
- To remove all the unsaved changes: git restore .
- To remove an unsaved file: git restore file_name

## Removing untracked changes:
- To remove untracked files: git clean -f
- To remove untracked folders: git clean -fd

## Adding modifications:
- To add all updated files: git add -u
- To add an updated file: git add file_name
- To commit the changes: git commit -m "message"
- To update the remote from local: git push

## Contributing to the official lineqGPR version:

- To push changes in the master : 
1. Push your changes in your working branch 
2. Make a pull request (PR) from Github to merge your branch into master    
3. Be sure the automatic tests and demos are validated before merging 

**Note:** This procedure must be followed in order to have a track of the changes in the master and allow peer code review

## Just in case:
- git reset --hard <commit_code>
- git push origin <your_branch_name> --force
- git branch --delete <branch_name>