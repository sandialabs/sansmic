{% if objtype == 'property' %}
:orphan:
{% endif %}

{{ objname | escape | underline}}

.. currentmodule:: {{ module }}

{% if objtype == 'property' %}
property
{% endif %}

.. auto{{ objtype }}:: {{ fullname | replace(module + ".", module + "::") }}
