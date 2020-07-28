.. _contributing:

Contributing
------------

We outline below the workflow we use for developing X-PSI.

Contact
~~~~~~~

You can interact with the X-PSI :ref:`team` in the public domain on GitHub.
We may also be contacted privately at *xpsi-team[at]googlegroups.com*.
Private development between releases is performed in a private repository on
BitBucket. We communicate internally in our Slack workspace.

To contribute to this project you can do so whilst being recognised as either
a community member or a team member. If you contribute to feature development
publicly on GitHub, you will be invited to be recognised as a team member in
the :ref:`acknowledgements`, and to join the mailing group and our Slack
workspace.

Community members can request to join our Slack workspace without being
recognised as a team member (or developing features). Community members who
are experienced users and are willing to donate some of their time to assisting
other X-PSI users (the potential benefits of engaging in such interaction are
self-explanatory) may also request to join the X-PSI team without developing
features.

Past developers are recognised as team members unless they no longer wish to be.
Aside from acknowledging your contribution, this is useful for feature
maintenance and a contact point for questions, if permission is granted.
Past team members will be acknowledged as such.

If you are a grant-holder and wish to base a major research project on the
development and/or application of X-PSI, we would like to invite you to discuss
your joining the X-PSI team in order to collaborate. If you find this idea
interesting then please contact Anna L. Watts, whose contact details are on the
:ref:`acknowledgements` page. You are entirely free to develop X-PSI only in the
public domain on GitHub, of course, but we can likely better organise and
execute this work by using a mixture of private and public platforms. If you
are applying X-PSI in research as an end-user, a private plaform for
communication is probably desirable, and we invite you to request access to our
Slack workspace.


.. _workflow:

Git workflow
~~~~~~~~~~~~

**GitHub repository**

X-PSI is an open-source software project and we aspire that all X-PSI versions
required for reproducibility of publications have release tags on GitHub.

End-users and community contributors can interact with the project freely on
GitHub. If you wish to contribute code, either in the form of a patch, a tweak,
or a feature, but are uncertain about the workflow, we now provide an advisory
workflow.

* Clone the ``xpsi`` repository to your local computer:

.. code-block:: bash

    git clone https://github.com/ThomasEdwardRiley/xpsi.git <path/to/xpsi>

* You will be on the ``master`` branch by default; this branch by default tracks
  the official (or central) ``origin/master`` branch, and moreover, ``master``
  is merely to act as an image of ``origin/master``, meaning that official
  upstream commits can be applied to ``master`` via fast-foward merge.

* Checkout a local branch for your work, which we'll assume is some patch, but
  could be a feature or otherwise:

.. code-block:: bash

    git checkout -b patch/fix_something

* Commit your work to ``patch/fix_something``:

.. code-block:: bash

    git commit -m 'patch something'

* Check to see if any there are any new upstream commits, which would mean
  that ``patch/fix_something`` and ``origin/master`` have diverged:

.. code-block:: bash

    git pull origin master:master

* If the branches have diverged, you can either rebase ``patch/fix_something``
  on ``master`` or merge ``master`` into ``patch/fix_something``, in either
  case resolving any conflicts:

.. code-block:: bash

    git rebase master <or> git merge master

* Note that if you have already pushed ``patch/fix_something`` to a remote
  repository you own (such as a fork of ``xpsi``; see below), and especially
  if this is accessible by others (e.g., via  submitted pull request), you
  should only consider merging ``master`` into ``patch/fix_something`` so the
  branch history is not rewritten.

* After integrating upstream changes, you might decide to continue to working
  on your branch or you could work on another branch for a different patch
  or feature, repeating the process of integrating upstream changes as
  appropriate, and as a requirement, to prepare for a pull request.

* Once you are ready to contribute your work to the ``xpsi`` repository,
  meaning that you have integrated any upsteam changes from ``xpsi``, you need
  a fork of the ``xpsi`` repository on the same hosting platform (GitHub). You
  can create a fork using the GitHub GUI.

* With the address of your ``fork`` you can add it as a remote to you local
  repository:

.. code-block:: bash

    git remote add fork https://github.com/<username>/xpsi.git

* Now push ``patch/fix_something`` to ``fork``, creating a remote branch
  ``fork/patch/fix_something`` that ``patch/fix_something`` tracks:

.. code-block:: bash

    git push -u fork

* Now you can submit a pull request, using the GitHub GUI, from
  ``fork/patch/fix_something`` to ``xpsi/master``.

* You can update the pull-request topic branch by pushing additional commits
  from ``patch/fix_something`` to ``fork/patch/fix_something``, which will
  update the pull request automatically:

.. code-block:: bash

    git push

* The pull request will then be reviewed and discussed, and will either be
  merged or merge will be pending, because we request additional commits on
  the pull-request topic branch. Your pull request may be declined in some
  instances because the work reproduces development work that has already been
  performed but not published; your pull request may also be ultimately
  declined if it contains changes or implementations that we do not support and
  which cannot for some reason be separated from changes we do support. Your
  intellectual contribution to the project will be gratefull acknowledged in
  the :ref:`acknowledgements` and/or in the project :ref:`history` if this
  interaction leads to some form of merged development/implementation by
  another community member even if your pull request is ultimately declined.

* If you co-authored a pull request with one or more collaborators, you can
  acknowledge them using the GitHub pull-request GUI as you would for a single
  commit. When a pull request is accepted, it is typically going to be via a
  merge-squash unless the history is clean or work will continue to be
  commited on the topic branch after the merge (where applicable). In this
  case it is the responsiblity of the X-PSI team member executing the merge
  to replicate the list of co-authors from the original pull request in the
  squash message.

* When a pull request is merged, conflicts will either need to be resolved
  locally by you as suggested above, ending in a pull request update, or by an
  X-PSI team member locally and then merged with or without a pull request.

* If you are ready to start development on a distinct patch or feature that is
  not conditional on your open pull requests being merged, then you can apply
  the workflow above by checking out a new local development branch off of
  an up-to-date ``master``. If your work *is* conditional on your open pull
  requests, you are free to do so by commiting to the relevant topic branch
  or with some other branching scheme, but there is a risk that more work will
  be needed if the open pull request is not merged into the central repository
  or only a subset of proposed changes are merged or conflict resolution does
  not favour all of the changes you proposed. Of course, such work may remain
  useful in your own applications even if it is never all merged into the
  central repository.

* The above workflow also applies to remote branches other than ``master`` that
  might exist in the ``xpsi`` repository that you wish to contribute to, but
  this will be a less common pattern.

If you want to contribute a feature, you are welcome to communicate with us
either on GitHub via issues and pull-requests, or on a private platform
(see below).


**BitBucket repository**

Most feature development by the X-PSI team is conducted on private platforms
including a private development repository.


. [#]_

.. rubric:: Footnotes

.. [#] .


