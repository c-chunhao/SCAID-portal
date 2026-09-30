<template>
  <main class="login-note">
    <h1>No user login</h1>
    <p>
      The public SCAID portal is read-only. Browse datasets, gene maps, and
      downloads without an account. There is no public password form and no
      token login API.
    </p>
    <p>
      <router-link to="/Home">Home</router-link>
      ·
      <router-link to="/Dataset">Dataset</router-link>
    </p>
  </main>
</template>

<script setup>
defineOptions({ name: 'LoginNotice' })
</script>

<style scoped>
.login-note {
  max-width: 40rem;
  margin: 4rem auto;
  padding: 0 1.25rem;
  color: #2c241c;
  font-family: Georgia, "Times New Roman", serif;
}
h1 {
  font-size: 1.75rem;
  margin-bottom: 0.75rem;
}
p {
  line-height: 1.5;
}
a {
  color: #9b5f47;
}
</style>
