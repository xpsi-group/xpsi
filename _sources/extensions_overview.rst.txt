.. _extensions:

Overview of extension modules
=============================

Users and developers are free to implemented functionality C extensions
to Python in order to increase speed. In the current framework some
functionality that a user is required to customize for model specification
is necessarily implemented via extension modules, whilst custom functionality
can optionally be written in extension modules.

Not all extension modules are documented yet. A subset of extension modules are
considered as reserved for internal usage, being accessed indirectly by
users through Python-exposed classes in the modeling language. Wrapping these
extension modules conceals lower-level handling of operations that the
user need not worry about: this part of the framework is less flexible because
the operations dominate likelihood function evaluation cost, and it is more
difficult to develop the source code whilst maintaining performance. The
extensions documented here are readily transplanable with other custom code,
be it extension modules or otherwise as appropriate for defining a custom
model, and in most cases should not have a substantial effect on overall speed.
