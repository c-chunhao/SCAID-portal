<template>
  <section class="container not-found">
    <h1>Page not found</h1>
    <p>This address does not match a SCAID page. Use the navigation above, or return to the dataset browser.</p>
    <router-link to="/Dataset">Browse datasets</router-link>
  </section>
</template>

<script setup>
defineOptions({ name: 'NotFoundPage' })
</script>

<style scoped>
.not-found { padding-top: 3rem; padding-bottom: 4rem; }
h1 { margin-bottom: 1rem; font-size: 2rem; }
p { max-width: 70ch; line-height: 1.6; }
a { display: inline-flex; align-items: center; min-height: 44px; color: var(--scaid-accent); text-underline-offset: .2em; }
</style>
