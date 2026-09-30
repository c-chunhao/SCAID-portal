<template>
  <div class="layout-container">
    <a class="skip-link" href="#main-content">Skip to main content</a>
    <LayoutHeader class="header" />
    <main id="main-content" class="main-content" tabindex="-1">
      <router-view />
    </main>
    <LayoutFooter class="footer" />
  </div>
</template>

<script setup>
defineOptions({ name: 'ScaidLayout' });
import LayoutHeader from './components/LayoutHeader.vue';
import LayoutFooter from './components/LayoutFooter.vue';

</script>

<style scoped>
.layout-container {
  display: flex;
  flex-direction: column;
  min-height: 100vh;
  min-height: 100dvh;
}

.header,
.footer {
  flex-shrink: 0;
}

.main-content {
  flex: 1;
  min-width: 0;
  scroll-margin-top: 88px;
}

.main-content:focus-visible {
  outline: 2px solid var(--scaid-accent, #8a503a);
  outline-offset: -2px;
}

.skip-link {
  position: fixed;
  top: 12px;
  left: 16px;
  z-index: 2100;
  padding: 12px 16px;
  border-radius: 6px;
  color: #fff;
  background: var(--scaid-accent, #8a503a);
  font-weight: 600;
  text-decoration: underline;
  text-underline-offset: 3px;
  transform: translateY(calc(-100% - 16px));
}

.skip-link:focus {
  transform: translateY(0);
  outline: 2px solid var(--scaid-ink, #2b2320);
  outline-offset: 3px;
}
</style>
