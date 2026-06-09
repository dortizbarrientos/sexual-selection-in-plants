# Git daily habits with one repo and two computers

Key rule to recall: **pull before you start, push before you stop.**
pull -> work -> commit -> push.  Never leave a machine without pushing.

## When I start working in a computer
    git pull --ff-only        # get whatever the other machine pushed
    # if it refuses ("divergent branches"):
    git pull --rebase         # replay your local commits on top (linear history)

## While I am working at a computer (commit in small, chunks, that you care about)
    git status                # what changed — run it constantly
    git add -A                # stage everything
    git commit -m "message"

## When I leave my computer
    git push
    # if rejected ("fetch first")
    git pull --rebase
    git push

## When I work with my other computer
    git stash                 # shelve uncommitted work
    git pull --ff-only
    git stash pop             # bring your work back on top

## Handy commands
    git status                # where am I / what changed
    git branch -vv            # branch + its remote tracking
    git diff                  # unstaged changes
    git diff --staged         # staged changes
    git log --oneline -10     # recent history
    git restore <file>        # discard unstaged edits to a file
    git restore --staged <f>  # unstage but keep edits
    git reset --soft HEAD~1   # undo last commit, keep changes
    git remote -v             # confirm the remote URL

## Two other things I should keep in mind when working
- Run `git status` before every `git add -A` so random zips never go up.
- And my favourite file of all times: Use `.gitignore` for stuff you want in
your computer but not in the cloud.  