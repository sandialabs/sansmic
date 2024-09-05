{{ objname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

{% if (methods or attributes) and not objname.endswith('Iterator') %}
   .. rubric:: Summary
   .. autosummary::

{% block methods %}
{% if methods %}
   {% for item in methods %}
      {%- if (not item.startswith('_') or item in ['__call__'] ) and item not in inherited_members %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}
   {% for item in attributes %}
      {%- if not item.startswith('_') and item not in inherited_members  %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
{% endif %}
{% endblock %}

   .. rubric:: Details

{% endif %}
